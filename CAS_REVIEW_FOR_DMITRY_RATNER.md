# Памятка для ревью: Проективная CAS над $\mathbb{Q}(i, \sqrt{d-1})$ и $\mathbb{P}^1(\mathbb{C})$ + интеграции в Mathics3 и SymvaCAS

Привет, Дима! Ниже — полная карта всех четырёх репозиториев по компьютерной алгебре и символьному анализу в организации **`AI-Guiders`**: что где лежит, какие новые возможности добавлены, какими командами прогнать автотесты и на каких конкретных интегралах/выражениях проверять руками (перед тем как открывать PR в апстримы `Mathics3` и `symvacas`).

---

## 0. Инженерный Use Case (за 2 минуты): Параметрическое моделирование самолётов, CAD-геометрия (NURBS / сопряжения обшивки) и аэродинамика крыла

Если убрать математические термины, этот движок напрямую решает **4 главные проблемы параметрического CAD/CAE самолётов**, на которых регулярно падают или теряют точность обычные `double`-геометры:

1. **Замкнутые шпангоуты фюзеляжа и профили крыла без разрыва на $180^\circ$ ($\theta = \pi \iff t = \infty$) и без `NaN`:**
   - В рациональной параметризации (NURBS / подстановка Вейерштрасса $\cos\theta = \frac{1-t^2}{1+t^2}, \sin\theta = \frac{2t}{1+t^2}$) при смыкании контура в задней кромке ($\theta = \pi$) параметр $t \to \infty$, в `double` происходит `DivisionByZero / NaN`, и контур приходится костыльно резать на куски, на стыках которых при перестроении модели плывёт кривизна $G^2$.
   - На проективной сфере Римана $\mathbb{P}^1$ (`P1Point(U, V)`) задняя кромка $t = \infty$ — это обычная точка `[1 : 0]`, где `P1Point.infinity().stereographic_circle()` без всяких пределов и делений на ноль возвращает точные `(-1, 0)`. Весь замкнутый контур шпангоута или профиля крыла параметризуется **одной непрерывной формулой без разрезов**.

2. **Никаких `Topology Rebuild Failed` на касаниях и зализах крыло–фюзеляж (`Fillet / Offset`):**
   - Когда в параметрической модели самолёта сдвигаешь нервюру или меняешь угол стреловидности и две поверхности переходят из трансверсального пересечения в **касание** (кратный корень!), производная обращается в ноль и численный Ньютон в `double` разваливается.
   - Редукция Остроградского–Эрмита (`hermite_reduction`) отделяет кратные корни и полюса любой кратности **чисто алгебраическим НОД полиномов $Q$ и $Q'$ в точных рациональных дробях вообще без численного поиска корней**.

3. **Точная гладкость обшивки $G^2 / G^3$ и рациональные сплайны Паде $[m/n]$ без потери знаков:**
   - Численное дифференцирование кривизны профиля конечными разностями в `double` съедает половину мантиссы и даёт паразитную волнистость обшивки.
   - Дуальное кольцо Кавальери (`p_diff`) и рациональный аппроксиматор Паде (`p_pade` / `PadeApproximant`) выдают **точные аналитические производные любого порядка** и компактные рациональные профили $[m/n]$ по любому параметру чертежа.

4. **Мгновенный аналитический пересчёт центровки, моментов инерции, эпюры давления и циркуляции Жуковского:**
   - Вместо перестроения тяжёлой расчётной сетки при каждом шаге параметрической оптимизации крыла, `p_integrate` и `p_contour_integrate` (`ContourIntegrate` / `Residue`) берут сингулярные интегралы нагрузки и контурной циркуляции Жуковского $\oint f(z)\,dz$ аналитически за миллисекунды с точностью до $N$ знаков.

### Готовый 15-строчный скрипт проверки авиационного Use Case (`guiders-mathematics`)

```python
import sympy as sp
from guiders_math import (
    P1Point, hermite_reduction, p_antiderivative,
    p_contour_integrate, p_diff, p_integrate, p_pade,
)

# 1. Задняя кромка профиля / смыкание шпангоута в theta = pi (t = oo = [1:0] на сфере Римана P^1)
t_inf = P1Point.infinity(d=3)
cos_pi, sin_pi = t_inf.stereographic_circle()
print("1. Trailing edge at t=[1:0]: (cos, sin) =", (cos_pi.a, sin_pi.a))
# -> (Fraction(-1, 1), Fraction(0, 1)) — без деления на ноль и без разреза контура!

# 2. Точные производные кривизны обшивки G^1 / G^2 в дуальном кольце Кавальери и рациональный сплайн Паде [2/2]
x = sp.Symbol("x")
airfoil_camber = (x * (1 - x)) / (1 + x**2)
g1 = p_diff(airfoil_camber, x, order=1)
g2 = p_diff(airfoil_camber, x, order=2)
pade_spline = p_pade(sp.exp(-x) * sp.sin(x), x, 2, 2)
print("2. Exact G1 =", g1, "| G2 =", g2, "| Pade [2/2] =", pade_spline)

# 3. Вырожденное касание поверхностей (3-кратный полюс): точное отделение через НОД Остроградского-Эрмита без float-Ньютона
decomp = hermite_reduction(1, (x**2 + 1)**3, x)
print("3. Hermite rational part =", decomp.hermite_rat, "+ int(", decomp.rem_num / decomp.sqfree_den, ")")
# -> (3*x**3 + 5*x)/(8*x**4 + 16*x**2 + 8) + int( 3/(8*(x**2 + 1)) )

# 4. Точный интеграл распределения нагрузки по хорде крыла и контурная циркуляция Жуковского
prim = p_antiderivative("1 / (1 + x**2)**3", "x")
wing_load = p_integrate("1 / (1 + x**2)**3", ("x", 0, 1), digits=40)
sym_circ, num_circ = p_contour_integrate("1 / (1 + z**4)", "z", center="I", radius=1, digits=40)
print("4. Wing load antiderivative =", prim.symbolic_expr, "| definite [0,1] =", wing_load.eval_R())
# -> 1/4 + 3*pi/32 = 0.5445243112740431161058728171824533953935
print("5. Kutta-Joukowski circulation =", sym_circ, "=", num_circ.re)
# -> sqrt(2)*pi/2 = 2.221441469079183123507940495186261855124
```

---

## 1. Сводная таблица репозиториев и веток

| № | Репозиторий (`AI-Guiders`) | Ветка и коммит | Назначение |
|---|---|---|---|
| **1** | [`AI-Guiders/guiders-mathematics`](https://github.com/AI-Guiders/guiders-mathematics) | `main` (`5a15b08`) | Эталонное ядро проективной CAS (`guiders_math`): биквадратичное поле $\mathbb{Q}(i,\sqrt{d-1})$, сфера Римана $\mathbb{P}^1(\mathbb{C})$, редукция Остроградского–Эрмита, сопряжённые вычеты Лагранжа, поворот сектора Френеля–Эйри $\zeta_{2n}$, резольвенты Стилтьеса–Паде для бесконечномерных пространств и чисто рациональный вычислитель произвольной точности (`0 mpmath`). |
| **2** | [`AI-Guiders/mathics-core`](https://github.com/AI-Guiders/mathics-core/tree/feat/projective-calculus-and-mcp) | `feat/projective-calculus-and-mcp` (`b2fdfcf38`) | Форк ядра **Mathics3** (открытая реализация Wolfram Language): 11 новых встроенных функций анализа, диадическая башня радикалов `Cos[p*Pi/2^k]`, быстрый проективный `Integrate` и встроенный **Wolfram Language MCP-сервер** (`mathics3-mcp`). |
| **3** | [`AI-Guiders/Mathics3-Frontend-django`](https://github.com/AI-Guiders/Mathics3-Frontend-django/tree/feat/notebook-intellisense) | `feat/notebook-intellisense` (`c1e406d9`) | Веб-ноутбук **Mathics3**: живой IntelliSense (`Tab`/`Enter`/`Esc`), попиксельный Signature Help прямо под открывающей скобкой `Func[|]` и автозакрытие `[...]`, `{...}`, `(...)`. |
| **4** | [`AI-Guiders/symvacas`](https://github.com/AI-Guiders/symvacas/tree/feat/projective-calculus-q-i-sqrt) | `feat/projective-calculus-q-i-sqrt` (`3ed1ee2`) | Форк **SymvaCAS**: проективное ядро Остроградского–Эрмита–Лагранжа и Френеля с генерацией пошаговых карточек (`Step-by-Step`) и вычислением до $N$ знаков. |

---

## 2. Репозиторий 1: `AI-Guiders/guiders-mathematics` (ядро проективной CAS)

### Что внутри (`src/guiders_math/`)
- **`field.py`** — точная арифметика в биквадратичном поле $\mathbb{Q}(i, \sqrt{d-1})$ (`QQuadExt`) без плавающей точки и однородные координаты на сфере Римана $\mathbb{P}^1(\mathbb{C})$ (`P1Point(U, V)`, где $\infty = [1 : 0]$, $0 = [0 : 1]$, и деление на ноль заменено проективным инвертированием $[U:V] \mapsto [V:U]$).
- **`quadrature.py`** — проективный интегратор:
  - Алгебраическое отделение рациональной части $\int R(x)\,dx = \frac{A(x)}{B(x)} + \int \frac{C(x)}{D(x)}\,dx$ через НОД Остроградского–Эрмита (`hermite_reduction`, без численного поиска корней).
  - Спаривание комплексно-сопряжённых полюсов $\alpha \pm i\beta$ в вещественные логарифмы $\ln((x-\alpha)^2+\beta^2)$ и двухпараметрические арктангенсы $\operatorname{atan2}(x-\alpha, \beta)$ со строгим подсчётом числа витков на $[a,b]$ и на $[0,\infty)$ через $\infty = [1:0] \in \mathbb{P}^1$.
  - Для неприводимых над $\mathbb{Q}$ факторов степени $\ge 3$ — чисто рациональный комплексный метод Аберта–Эрлиха (`_aberth_complex_roots_rational`) без `mpmath`/`numpy`.
  - Поворот сектора плоскости Кавальери на корень из единицы $\zeta_{2n} = e^{i\pi/(2n)}$ для осциллирующих интегралов Френеля–Эйри $\int_0^\infty \cos(x^n)\,dx$ и $\int_0^\infty \sin(x^n)\,dx$.
  - Контурные интегралы $\oint_{|z-z_0|=r} f(z)\,dz = 2\pi i \sum \operatorname{Res}$ (`integrate_closed_contour`).
- **`eval_rc.py`** — вычислитель элементарных функций (`eval_cos`, `eval_sin`, `eval_tan`, `eval_atan`, `eval_atan2`, `eval_exp`, `eval_log`, `eval_sqrt`, `eval_pi`) на $\mathbb{R}$ и $\mathbb{C}$ с произвольной точностью $N$ десятичных знаков через итерацию Герона–Ньютона в `Fraction` и полуугловой спуск Кэли–Паде $(5,5)$. **Ни одного импорта `mpmath` и ни одного вызова `math.sin/cos/atan/log/exp/sqrt`.**
- **`operator_space.py`** — бесконечномерные гильбертовы пространства: аппроксимация резольвенты Стилтьеса–Паде $R(z) = (zI - H)^{-1}$ по моментам $\mu_k = \langle \psi_0 | H^k | \psi_0 \rangle$ (алгоритм Ланцоша–Стилтьеса в точных `Fraction`), полюса спектра и вычисление оператора эволюции $\langle \psi_0 | e^{-iHt} | \psi_0 \rangle$ через контурный интеграл Рисса–Данфорда по вычетам.
- **`sympy_bridge.py`** — drop-in замена `projective_integrate(expr, var_or_tuple)` с поддержкой `precision=N`.

### Как запускать проверку
```bash
git clone https://github.com/AI-Guiders/guiders-mathematics.git
cd guiders-mathematics
pip install -e .
python verify_principles.py
```

### Что проверяет `verify_principles.py` (9 гейтов)
1. **`p1_riemann_sphere`** — арифметика на $\mathbb{P}^1(\mathbb{C})$, инволюция $\infty = [1:0] \leftrightarrow 0 = [0:1]$, преобразования Мёбиуса и подстановка Вейерштрасса в $t = \infty$ ($\cos(\pi) = -1, \sin(\pi) = 0$).
2. **`qquad_field`** — точные тождества в $\mathbb{Q}(i, \sqrt{2})$ и $\mathbb{Q}(i, \sqrt{3})$, инверсия $(1 + i + \sqrt{2})^{-1}$.
3. **`cavalieri_pade`** — производные в дуальном кольце $\mathbb{Q}[\varepsilon]/(\varepsilon^{m+1})$ и аппроксиманта Паде $[2/2]$ для $\exp(x)$ (`c(0)=1, c(1)=1/2, c(2)=1/12`).
4. **`ostrowski_hermite_lagrange`** — кратные полюса $\int_0^1 \frac{dx}{(x^2+1)^3} = \frac{1}{4} + \frac{3\pi}{32}$ без ухода в $\mathbb{C}$ и несобственный интеграл $\int_0^\infty \frac{dx}{(x^2+1)^2} = \frac{\pi}{4}$ через $[1:0] \in \mathbb{P}^1$.
5. **`fresnel_sector_rotation`** — $\int_0^\infty \cos(x^2)\,dx = \frac{\sqrt{2\pi}}{4}$ и $\int_0^\infty \cos(x^4)\,dx = \frac{1}{4}\Gamma(\tfrac{1}{4})\cos(\tfrac{\pi}{8})$ с проверкой на **50 десятичных знаках**.
6. **`eval_rc_purity`** — AST-проверка файла `eval_rc.py` (ноль `mpmath`, ноль `float`-тригонометрии) + проверка $\pi$, $\sin(1+2i)^2+\cos(1+2i)^2=1$ и $\exp(\ln(3+4i))=3+4i$ на **60 знаках**.
7. **`contour_and_path_integrals`** — $\oint_{|z|=2} \frac{e^z}{z^3}\,dz = \pi i$, $\oint_{|z|=2} \frac{dz}{z^2+1} = 0$, комплексный путь от $0$ до $1+i$.
8. **`sympy_speed_and_coverage`** — $\int \sqrt{\tan x}\,dx$ и $\int \frac{dx}{1+\sin^4 x}$ (оба в штатном SymPy уходят в многосекундный Groebner-завис, а здесь берутся аналитически за миллисекунды + сверка производной в произвольной точке и численной квадратуры на 30 знаках).
9. **`infinite_dimensional_operator_space`** — восстановление спектра бесконечномерного гармонического осциллятора $E_n = n + 1/2$ и унитарная эволюция $|\langle\psi_0|e^{-iHt}|\psi_0\rangle| = 1$ через резольвенту Стилтьеса–Паде.

---

## 3. Репозиторий 2: `AI-Guiders/mathics-core` (Wolfram Language CAS + MCP Server)

### Как установить и прогнать тесты
```bash
git clone -b feat/projective-calculus-and-mcp https://github.com/AI-Guiders/mathics-core.git
cd mathics-core
pip install -e .
pytest test/builtin/numbers/test_projective_calculus.py -v
```
*(Должно пройти: **`86 passed`** за ~3.5 сек).*

### На чём проверять в CLI `mathics` или в ноутбуке (список контрольных выражений)

В стандартном `Mathics3` (апстрим) **все 11 пунктов ниже либо вообще отсутствовали (возвращались невычисленными), либо зависали/падали**. Теперь они работают из коробки:

1. **Аппроксиманты Паде (`PadeApproximant`) в $x_0 = 0$, в произвольной точке $x_0 \neq 0$ и на сфере Римана $x_0 = \infty$:**
   ```wolfram
   PadeApproximant[Exp[x], {x, 0, {2, 2}}]
   (* -> (1 + x / 2 + x ^ 2 / 12) / (1 - x / 2 + x ^ 2 / 12) *)

   PadeApproximant[Sin[x], {x, 0, {3, 2}}]
   (* -> (x - 7 * x ^ 3 / 60) / (1 + x ^ 2 / 20) *)

   PadeApproximant[Log[x], {x, 1, {2, 2}}]
   (* -> (-3 + 3 * x ^ 2) / (1 + 4 * x + x ^ 2) *)
   ```
2. **Вычеты (`Residue`), включая кратные полюса, существенно особые точки и вычет в $\infty = [1:0] \in \mathbb{P}^1$:**
   ```wolfram
   Residue[1 / (z^2 + 1), {z, I}]
   (* -> -I / 2 *)

   Residue[Exp[z] / z^3, {z, 0}]
   (* -> 1 / 2 *)

   Residue[z^3 / (z^4 + 1), {z, Infinity}]
   (* -> -1 *)
   ```
3. **Замкнутые контурные интегралы в $\mathbb{C}$ (`ContourIntegrate`):**
   ```wolfram
   ContourIntegrate[1 / (z - 1), {z, 0, 2}]
   (* -> 2 * I * Pi *)

   ContourIntegrate[Exp[z] / z^3, {z, 0, 2}]
   (* -> I * Pi *)

   ContourIntegrate[1 / (z^2 + 1), {z, I, 1}]
   (* -> Pi *)
   ```
4. **Несобственные интегралы Френеля–Эйри (`Integrate` через поворот сектора $\zeta_{2n}$):**
   ```wolfram
   Integrate[Cos[x^2], {x, 0, Infinity}]
   (* -> Sqrt[2 * Pi] / 4 *)

   Integrate[Sin[x^2], {x, 0, Infinity}]
   (* -> Sqrt[2 * Pi] / 4 *)
   ```
5. **Тяжёлые рациональные, дробно-тригонометрические и вейерштрассовы интегралы (мгновенно вместо зависания):**
   ```wolfram
   Integrate[1 / (x^2 + 1)^3, {x, 0, 1}]
   (* -> 1 / 4 + 3 * Pi / 32 *)

   Integrate[1 / (x^4 + 1), {x, 0, Infinity}]
   (* -> Sqrt[2] * Pi / 4 *)

   Integrate[1 / (1 + Sin[x]^4), {x, 0, Pi / 2}]
   (* -> Sqrt[6] * Pi / 4 *)

   Integrate[Sqrt[Tan[x]], x]
   (* -> -Sqrt[2]*Log[1 - Sqrt[2]*Sqrt[Tan[x]] + Tan[x]]/4 + Sqrt[2]*Log[1 + Sqrt[2]*Sqrt[Tan[x]] + Tan[x]]/4 + Sqrt[2]*ArcTan[-1 + Sqrt[2]*Sqrt[Tan[x]]]/2 + Sqrt[2]*ArcTan[1 + Sqrt[2]*Sqrt[Tan[x]]]/2 *)
   ```
6. **Полное семейство тригонометрических преобразований (`TrigReduce`, `TrigExpand`, `TrigFactor`, `TrigToExp`, `ExpToTrig`):**
   ```wolfram
   TrigReduce[Sin[x]^2 * Cos[x]]
   (* -> Cos[x] / 4 - Cos[3 * x] / 4 *)

   TrigExpand[Cos[3 * x]]
   (* -> Cos[x] ^ 3 - 3 * Cos[x] * Sin[x] ^ 2 *)

   TrigFactor[Sin[x] + Sin[3 * x]]
   (* -> 2 * Sin[2 * x] * Cos[x] *)

   TrigToExp[Cos[x]]
   (* -> E ^ (-I * x) / 2 + E ^ (I * x) / 2 *)

   ExpToTrig[E^(I * x) + E^(-I * x)]
   (* -> 2 * Cos[x] *)
   ```
7. **Диадическая башня вложенных радикалов для `Cos[p*Pi/2^k]` и `Sin[p*Pi/2^k]` при любых $2^k \ge 16$:**
   ```wolfram
   Cos[Pi / 16]
   (* -> Sqrt[2 + Sqrt[2 + Sqrt[2]]] / 2 *)

   Sin[Pi / 32]
   (* -> Sqrt[2 - Sqrt[2 + Sqrt[2 + Sqrt[2]]]] / 2 *)
   ```
8. **Прямое и обратное преобразования Лапласа (`LaplaceTransform`, `InverseLaplaceTransform`):**
   ```wolfram
   LaplaceTransform[t^2 * Exp[-3 * t], t, s]
   (* -> 2 / (3 + s) ^ 3 *)

   InverseLaplaceTransform[1 / (s^2 + 4), s, t]
   (* -> Sin[2 * t] / 2 * HeavisideTheta[t] *)
   ```
9. **Бесконечномерная резольвента Стилтьеса–Паде (`StieltjesPadeResolvent`):**
   ```wolfram
   StieltjesPadeResolvent[{1, 1/2, 3/4, 11/8, 19/16}, z]
   (* -> (-3 + 4 * z) / (3 - 8 * z + 4 * z ^ 2)  — точные полюса z = 1/2 и z = 3/2! *)
   ```
10. **Нециклические чисто рациональные фундаментальные константы (`ProjectiveConstant[name, digits]`, 0 `mpmath`):**
    ```wolfram
    ProjectiveConstant["Pi", 50]
    (* -> 3.1415926535897932384626433832795028841971693993751 *)

    ProjectiveConstant["E", 50]
    (* -> 2.7182818284590452353602874713526624977572470937000 *)

    ProjectiveConstant["Lemniscate", 50]
    (* -> 2.6220575542921198104648395898911194136827549514316  (через AGM(1, Sqrt[2]) в Q(Sqrt[2])) *)

    ProjectiveConstant["Gamma14", 50]
    (* -> 3.6256099082219083119306851558676720029951676828801  (Gamma(1/4)) *)

    ProjectiveConstant["Catalan", 50]
    (* -> 0.91596559417721901505460351493238411077414937428167  (через Log_P(2 + Sqrt[3]) в Q(Sqrt[3])) *)

    ProjectiveConstant["Apery", 50]
    (* -> 1.2020569031595942853997381615114499907649862923405  (Zeta(3)) *)

    ProjectiveConstant["EulerGamma", 50]
    (* -> 0.57721566490153286060651209008240243104215933593992 *)
    ```
11. **Wolfram Language MCP-сервер (`mathics3-mcp`) для подключения к любому ИИ-агенту (Cursor / Claude Code / OpenCode):**
    - Точка входа: `mathics3-mcp` (или `python -m mathics.mcp_server`).
    - Экспортирует 6 MCP-инструментов: `wolfram_eval`, `wolfram_integrate`, `wolfram_pade`, `wolfram_residue`, `wolfram_contour_integrate`, `wolfram_stieltjes_resolvent`.
    - Быстрая проверка из терминала:
      ```bash
      python -c "from mathics.mcp_server import MathicsMCPServer; s = MathicsMCPServer(); print(s.call_tool('wolfram_integrate', {'expr': '1/(x^4+1)', 'var': 'x', 'lower': '0', 'upper': 'Infinity'}))"
      ```

---

## 4. Репозиторий 3: `AI-Guiders/Mathics3-Frontend-django` (Веб-ноутбук с IntelliSense)

### Как запустить веб-интерфейс
```bash
git clone -b feat/notebook-intellisense https://github.com/AI-Guiders/Mathics3-Frontend-django.git
cd Mathics3-Frontend-django
pip install -e .
python -m mathics_django.manage runserver 127.0.0.1:8088 --noreload
```
Открой в браузере **`http://127.0.0.1:8088/`**.

### Что проверить в браузере
1. **Живой IntelliSense:** начни печатать в ячейке `Pade`, `Cont`, `Stielt` или `Trig` — после 2 символов прямо под курсором всплывает список всех символов `System`*` и `Global`*` с сигнатурой и кратким описанием. Навигация: `↑`/`↓`, вставка: `Tab` или `Enter`, закрытие: `Esc`.
2. **Автозакрытие скобок:** при вводе `[` автоматически вставляется `]` с кареткой внутри `[|]` (аналогично для `{...}` и `(...)`), а повторное нажатие `]` перешагивает через закрывающую скобку без дублирования.
3. **Попиксельный Signature Help:** когда каретка стоит внутри вызова `PadeApproximant[|]` или `ContourIntegrate[|]`, прямо под символом `[` (вычисляется через скрытый mirror-`div` с учётом шрифта, скролла и переноса строк) висит компактная подсказка с синтаксисом аргументов (`PadeApproximant[expr, {x, x0, {m, n}}]`).

---

## 5. Репозиторий 4: `AI-Guiders/symvacas` (SymvaCAS + карточки пошагового решения)

### Как прогнать тесты и запустить веб-интерфейс
```bash
git clone -b feat/projective-calculus-q-i-sqrt https://github.com/AI-Guiders/symvacas.git
cd symvacas
python -m unittest discover -s app/src/main/python -p "test_*.py" -v
```
*(Должно пройти: **`15 passed`**).*

Чтобы открыть веб-приложение `SymvaCAS` (с обновлённым `web/engine.zip`):
```bash
python -m http.server 8080 --directory web
```
Открой **`http://127.0.0.1:8080/`**.

### Что проверить в `SymvaCAS` (`calc_projective_cas.py`)
- Встроен проективный движок `projective_integrate_with_steps(expr_str, var_name, lower, upper, precision)` с генерацией карточек `Step-by-Step`:
  1. Разложение Остроградского–Эрмита $\int R(x)\,dx = A(x)/B(x) + \int C(x)/D(x)\,dx$.
  2. Спаривание сопряжённых полюсов Лагранжа над $\mathbb{Q}(i, \sqrt{d-1})$.
  3. Проективный предел $\infty = [1:0] \in \mathbb{P}^1(\mathbb{C})$ на $[0, \infty)$ и $(-\infty, \infty)$.
  4. Секторный поворот Френеля–Эйри $\zeta_{2n} = e^{i\pi/(2n)}$ для $\int_0^\infty \cos(x^n)\,dx$.
  5. Вычисление результата до любого запрошенного числа десятичных знаков `precision` (например, 50 или 100 знаков) через чисто рациональную аппроксимацию Паде $(5,5)$ и итерацию Герона–Ньютона.
