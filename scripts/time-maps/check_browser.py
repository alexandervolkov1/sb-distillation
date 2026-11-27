"""Inspect the generated Clay notebook in Edge and capture verification screenshots."""

import argparse
from pathlib import Path

from playwright.sync_api import sync_playwright


CAPTURE_VEGA = """
window.mapsQA = [];
let currentEmbed;
Object.defineProperty(window, 'vegaEmbed', {
  configurable: true,
  get: () => currentEmbed,
  set: fn => { currentEmbed = (...args) => fn(...args).then(result => {
    if (args[0].closest?.('.time-map')) {
      window.mapsQA.push({element: args[0], spec: args[1], view: result.view});
    }
    return result;
  }); }
});
"""

FITS = """() => {
  const cards = [...document.querySelectorAll('.time-map')];
  return cards.length === 8 && cards.every(card => {
    const chart = card.querySelector('svg.marks');
    return chart && chart.getBoundingClientRect().right <= card.getBoundingClientRect().right + 0.5;
  });
}"""

VERIFY_SPECS = """() => {
  const reports = window.mapsQA.map(({element, spec, view}) => {
    const card = element.closest('.time-map');
    const feed = Number(card.closest('.time-map-row').dataset.feed);
    const {encoding, data} = spec.layer[0];
    const levels = spec.layer.slice(1, 3).flatMap(layer => layer.data.values.map(row => row.level));
    const groups = view.scenegraph().root.items[0].items;
    const rects = groups.find(item => item.name === 'layer_0_marks').items;
    const svg = element.querySelector('svg.marks');
    const titles = [...svg.querySelectorAll('.role-axis-title')].map(label => label.textContent);
    const ticks = [...svg.querySelectorAll('.role-legend-label text')].map(label => label.textContent);
    const width = svg.getBoundingClientRect().width;
    return {
      model: card.dataset.model, feed, count: rects.length, color: encoding.color.scale,
      width, titles, ticks, x: encoding.x.scale.domain, y: encoding.y.scale.domain,
      levels: [...new Set(levels)].sort((a, b) => a - b),
      min: Math.min(...data.values.map(d => d.days)),
      max: Math.max(...data.values.map(d => d.days)),
      text: spec.layer.at(-1).data.values.length,
      cells: feed === 0 && data.values.map(({model, ...cell}, index) => [cell, rects[index].fill])
    };
  });
  const zero = reports.filter(r => r.feed === 0);
  return {reports: reports.map(({cells, ...item}) => item),
          zeroEqual: JSON.stringify(zero[0].cells) === JSON.stringify(zero[1].cells)};
}"""


def check_page(file, output, channel):
    faults = []
    with sync_playwright() as driver, driver.chromium.launch(channel=channel) as browser:
        page = browser.new_page(viewport={"width": 1440, "height": 960}, ignore_https_errors=True)
        page.on("console", lambda message: faults.append(message.text)
                if message.type in ("error", "warning") else None)
        page.on("pageerror", lambda error: faults.append(str(error)))
        page.add_init_script(CAPTURE_VEGA)
        page.goto(file.resolve().as_uri(), wait_until="load", timeout=90000)
        page.wait_for_function("window.mapsQA.length === 8")
        page.wait_for_function(FITS)
        page.evaluate("document.fonts.ready")
        assert page.evaluate("""() => [...document.fonts].some(f => f.family === 'Old Standard'
                                                                 && f.status === 'loaded')""")
        inspection = page.evaluate(VERIFY_SPECS)
        reports = inspection["reports"]
        assert inspection["zeroEqual"]
        assert {(r["feed"], r["model"]) for r in reports} == {
            (feed, model) for feed in (0, 100, 200, 300) for model in ("cpim", "dgm")
        }
        domain = [min(r["min"] for r in reports), max(r["max"] for r in reports)]
        for report in reports:
            assert report["count"] == 2501
            assert report["color"]["type"] == "log" and report["color"]["domain"] == domain
            assert abs(report["width"] - reports[0]["width"]) <= 1
            assert report["x"] == [6, 12] and report["y"] == [680, 720]
            assert {"Внутренний диаметр, мм", "Температура источника, °C"} <= set(report["titles"])
            assert {"1", "2", "4", "8", "16", "32"} <= set(report["ticks"])
            expected = [v for v in [0.25, 0.5, 0.75, *range(1, int(domain[1]) + 1)]
                        if report["min"] < v < report["max"]]
            assert report["levels"] == expected and report["text"] > 0
        assert page.locator("pre:visible").count() == 0
        assert page.evaluate("document.documentElement.scrollWidth === innerWidth")
        output.mkdir(parents=True, exist_ok=True)
        page.locator(".distillation-maps").screenshot(path=output / "maps-overview.png")
        page.screenshot(path=output / "notebook-full.png", full_page=True)
        pair = page.locator(".time-map-row[data-feed='200']")
        pair.scroll_into_view_if_needed()
        target = page.evaluate("""() => {
            const {element, view} = window.mapsQA.find(({element}) => {
              const card = element.closest('.time-map');
              return card.dataset.model === 'dgm' && card.closest('.time-map-row').dataset.feed === '200';
            });
            const rects = view.scenegraph().root.items[0].items.find(i => i.name === 'layer_0_marks').items;
            const index = rects.findIndex(({datum}) => datum['temp-c'] === 700 && datum['diameter-mm'] === 10);
            const rect = element.querySelector('.layer_0_marks').children[index].getBoundingClientRect();
            return [rect.left + rect.width / 2, rect.top + rect.height / 2];
        }""")
        page.mouse.move(*target)
        tooltip = page.locator("#vg-tooltip-element.visible")
        tooltip.wait_for()
        tooltip_text = tooltip.inner_text()
        for expected in ("DGM", "200", "700", "10.0", "8.307", "60.19", "2.508"):
            assert expected in tooltip_text, tooltip_text
        page.screenshot(path=output / "maps-tooltip.png")
        page.set_viewport_size({"width": 390, "height": 900})
        page.wait_for_function(FITS)
        assert page.evaluate("document.documentElement.scrollWidth === innerWidth")
        stacked = pair.locator(".time-map").evaluate_all(
            "cards => cards[1].getBoundingClientRect().top >= cards[0].getBoundingClientRect().bottom"
        )
        assert stacked
        pair.screenshot(path=output / "maps-mobile.png")
        assert not faults, faults
        print(f"Eight charts and 20008 cells, common domain {domain}")
        print("Desktop and mobile layout, tooltip, contours, labels, font and hidden code verified.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("file", type=Path)
    parser.add_argument("screenshots", type=Path)
    parser.add_argument("--channel", default="msedge")
    options = parser.parse_args()
    check_page(options.file, options.screenshots, options.channel)
