(ns sb-distillation.diagrams.apparatus)

(defn boundary-conditions-svg
  "Black-and-white diagram of the CPIM apparatus
  and pressure boundary conditions."
  []
  (str
   "<svg xmlns='http://www.w3.org/2000/svg'
         viewBox='0 0 720 465'
         width='100%'
         font-family='Old Standard, serif'
         font-size='20'
         fill='black'>

      <defs>
        <marker id='arrow' markerWidth='8' markerHeight='8'
                refX='7' refY='3' orient='auto-start-reverse'>
          <path d='M0 0 L7 3 L0 6 Z' fill='black'/>
        </marker>

        <pattern id='hatch' width='6' height='6'
                 patternUnits='userSpaceOnUse'>
          <path d='M0 6 L6 0' stroke='black' stroke-width='0.7'/>
        </pattern>
      </defs>

      <rect width='720' height='465' fill='white'/>

      <!-- Ampoules and connecting tube -->
      <g fill='none' stroke='black' stroke-width='2'>
        <path d='M210 155 V125 Q210 90 175 90
                 H60 Q25 90 25 125 V265
                 Q25 300 60 300 H175 Q210 300 210 265 V205'/>

        <path d='M510 155 V125 Q510 90 545 90
                 H660 Q695 90 695 125 V265
                 Q695 300 660 300 H545 Q510 300 510 265 V205'/>

        <path d='M210 155 H510 M210 205 H510'/>
      </g>

      <!-- Liquid antimony -->
      <path d='M26 250 Q117 240 209 250 V265
               Q209 299 175 299 H60 Q26 299 26 265 Z'
            fill='url(#hatch)' stroke='black'/>

      <!-- Deposited crystals -->
      <g fill='white' stroke='black' stroke-width='1.5'>
        <path d='M545 298 L535 278 L545 261 L555 278 Z'/>
        <path d='M580 298 L569 267 L580 248 L591 267 Z'/>
        <path d='M615 298 L605 280 L615 264 L625 280 Z'/>
        <path d='M650 298 L640 270 L650 252 L660 270 Z'/>
      </g>

      <!-- Gas inlet and outlet -->
      <g stroke='black' stroke-width='1.5'
         fill='none' marker-end='url(#arrow)'>
        <path d='M117 25 V90'/>
        <path d='M602 90 V45 H690'/>
      </g>

      <text x='132' y='42'>H₂ от РРГ</text>
      <text x='650' y='30' text-anchor='middle'>К ловушке</text>

      <!-- Species fluxes -->
      <g stroke='black' stroke-width='1.5'
         marker-end='url(#arrow)'>
        <line x1='235' y1='170' x2='485' y2='170'/>
        <line x1='235' y1='190' x2='485' y2='190'/>
      </g>

      <text x='360' y='140' text-anchor='middle'>N₁, N₂</text>

      <!-- Apparatus labels -->
      <g text-anchor='middle'>
        <text x='117' y='120'>Источник</text>
        <text x='602' y='120'>Приёмник</text>
        <text x='117' y='225'>Sb (жидк.)</text>
        <text x='602' y='225'>Sb (тв.)</text>
      </g>

      <!-- Tube length and coordinate references -->
      <g stroke='black' stroke-width='1'
         stroke-dasharray='4 4'>
        <line x1='210' y1='215' x2='210' y2='350'/>
        <line x1='510' y1='215' x2='510' y2='350'/>
      </g>

      <line x1='210' y1='335' x2='510' y2='335'
            stroke='black'
            marker-start='url(#arrow)'
            marker-end='url(#arrow)'/>

      <g text-anchor='middle'>
        <text x='360' y='325'>L</text>
        <text x='210' y='365'>z = 0</text>
        <text x='510' y='365'>z = L</text>

        <text x='117' y='392'>p₁(0) = p₁,sat(Tₛ)</text>
        <text x='117' y='420'>p₂(0) = ?</text>

        <text x='602' y='392'>p₁(L) = 0</text>
        <text x='602' y='420'>p₂(L) = p_g + Q₂/C_trap</text>
      </g>

      <!-- Backward integration -->
      <line x1='495' y1='392' x2='225' y2='392'
            stroke='black' stroke-width='1.5'
            marker-end='url(#arrow)'/>

      <text x='360' y='420' text-anchor='middle'
            font-size='16'>Интегрирование</text>

    </svg>"))
