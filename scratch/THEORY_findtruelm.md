# Extraction des exposants de Floquet d'un spectre HSS : théorie, difficultés, solution

Document de travail accompagnant `findTruelmV3.m`, `findTrueFloquetV3.m` et la suite de tests `stress_v3.m` / `test_findTruelmV3.m`. Il pose le cadre mathématique, isole chaque difficulté numérique, et justifie l'algorithme retenu.

---

## 1. Cadre et objectif

On considère un système linéaire périodique (LTP) en temps continu

$$
\dot{x}(t) = A(t)\,x(t), \qquad A(t+T) = A(t), \qquad \omega = \frac{2\pi}{T},
$$

avec $A(t) \in \mathbb{R}^{n_x \times n_x}$. La stabilité et les modes propres d'un tel système ne se lisent pas sur les valeurs propres de $A(t)$ à un instant donné : ils sont gouvernés par la théorie de Floquet. L'objectif pratique est d'extraire les $n_x$ **exposants de Floquet** $\mu_k$ à partir du spectre d'une représentation harmonique tronquée (Hill / Harmonic State Space, HSS).

Ce problème apparaît dès qu'on analyse la stabilité d'un convertisseur de puissance en régime établi, d'une machine en régime déséquilibré, ou de toute dynamique linéarisée autour d'une trajectoire périodique.

---

## 2. Théorie de Floquet

Soit $\Phi(t)$ la matrice fondamentale de $\dot{x}=A(t)x$, avec $\Phi(0)=I$. Le théorème de Floquet établit la factorisation

$$
\Phi(t) = P(t)\,e^{R t}, \qquad P(t+T)=P(t), \quad P(0)=I,
$$

où $P$ est périodique et $R$ constante. La **matrice de monodromie** est

$$
M = \Phi(T) = e^{R T}.
$$

Ses valeurs propres $\rho_k$ sont les **multiplicateurs de Floquet**, et les **exposants de Floquet** $\mu_k$ sont définis par

$$
\rho_k = e^{\mu_k T}.
$$

La stabilité asymptotique équivaut à $|\rho_k|<1$ pour tout $k$, soit $\operatorname{Re}\mu_k<0$.

### 2.1 Les exposants vivent sur un cylindre

Le logarithme complexe n'est pas injectif : si $\rho_k = e^{\mu_k T}$, alors pour tout $m\in\mathbb{Z}$,

$$
e^{(\mu_k + j m\omega)T} = e^{\mu_k T}\,e^{j 2\pi m} = \rho_k .
$$

Les exposants ne sont donc définis qu'à un multiple entier de $j\omega$ près. Ils appartiennent au cylindre

$$
\mu_k \in \mathbb{C}\,/\,(j\omega\,\mathbb{Z}),
$$

c'est-à-dire un plan complexe dont l'axe imaginaire est replié périodiquement sur l'intervalle $[-\omega/2,\,\omega/2)$ (zone de Brillouin). **Ce caractère circulaire de l'axe imaginaire conditionne toute la suite.**

---

## 3. Représentation HSS et structure en peigne

On cherche les solutions de Floquet sous la forme $x(t) = e^{\mu t}\sum_{n} x_n\, e^{j n\omega t}$. En injectant dans $\dot{x}=A(t)x$ et en équilibrant chaque harmonique, on obtient un problème aux valeurs propres de dimension infinie sur les coefficients $\{x_n\}$. L'opérateur associé est un opérateur de Toeplitz par blocs (les harmoniques $A_{m-n}$ de $A(t)$) auquel s'ajoute le terme diagonal $-jn\omega I$.

La troncature aux harmoniques $|n|\le h$ donne la **matrice de Hill** $H_h$, de taille $(2h+1)\,n_x$. Ses valeurs propres approchent l'ensemble

$$
\big\{\, \mu_k + j m\omega \;:\; k=1,\dots,n_x,\; m=-h,\dots,h \,\big\}.
$$

Chaque exposant vrai $\mu_k$ engendre donc un **peigne** d'alias régulièrement espacés de $\omega$ sur l'axe imaginaire. Une fois l'axe imaginaire replié dans la zone de Brillouin, tout le peigne s'effondre sur un point unique $(\operatorname{Re}\mu_k,\ \operatorname{wrap}(\operatorname{Im}\mu_k))$.

**Le problème d'extraction se formule ainsi :** à partir des $(2h+1)\,n_x$ valeurs propres de $H_h$, identifier les $n_x$ exposants vrais, c'est-à-dire un représentant par peigne, en rejetant les alias et la pollution numérique.

---

## 4. Multiplicité, matrices défectives, points exceptionnels

L'extraction se complique lorsque plusieurs exposants coïncident. Les définitions suivantes fixent le vocabulaire.

### 4.1 Multiplicité algébrique et géométrique

Pour une valeur propre $\rho$ d'une matrice $M$ :

- la **multiplicité algébrique** $a(\rho)$ est sa multiplicité comme racine du polynôme caractéristique $\det(M-\rho I)$ ;
- la **multiplicité géométrique** $g(\rho) = \dim\ker(M-\rho I)$ est le nombre de vecteurs propres indépendants associés.

On a toujours $1 \le g(\rho) \le a(\rho)$.

### 4.2 Matrice défective

Une matrice est **défective** si $g(\rho) < a(\rho)$ pour au moins une valeur propre, c'est-à-dire si elle n'est pas diagonalisable. Sa forme de Jordan contient alors un bloc de taille $>1$,

$$
J_m(\rho) = \begin{pmatrix} \rho & 1 & & \\ & \rho & \ddots & \\ & & \ddots & 1 \\ & & & \rho \end{pmatrix} \in \mathbb{C}^{m\times m},
$$

de multiplicité algébrique $m$ mais géométrique $1$.

### 4.3 Conséquence numérique : éclatement en $\varepsilon^{1/m}$

Le point central pour notre problème est le **mauvais conditionnement spectral** des matrices défectives. Une perturbation de norme $\varepsilon$ appliquée à un bloc $J_m(\rho)$ déplace sa valeur propre $m$-fois dégénérée en $m$ valeurs réparties sur un cercle de rayon

$$
|\Delta\rho| \sim \varepsilon^{1/m}.
$$

En arithmétique flottante, $\varepsilon \approx 2\times 10^{-16}$, donc un mode défectif n'apparaît jamais comme $m$ valeurs propres rigoureusement égales : il **éclate** d'un montant

$$
m=2 \Rightarrow \sim 10^{-8}, \qquad m=3 \Rightarrow \sim 6\times 10^{-6}, \qquad m=5 \Rightarrow \sim 10^{-3}.
$$

Mesuré sur des blocs de Jordan conjugués par une similarité pleine, l'éclatement observé suit cette loi (de $5\times10^{-8}$ pour $m=2$ à $4\times10^{-3}$ pour $m=5$).

### 4.4 Point exceptionnel

Un **point exceptionnel** (EP) est une valeur du paramètre où deux exposants *et* leurs vecteurs propres coalescent simultanément, produisant un bloc de Jordan. En théorie de Floquet, les EP se situent génériquement sur les bords des langues de stabilité (équation de Mathieu / Hill). Au voisinage d'un EP, l'écart entre les deux exposants se comporte en $\sqrt{\delta}$ où $\delta$ mesure la distance au point : les modes sont proches mais distincts, puis fusionnent exactement à l'EP.

---

## 5. Les difficultés, une par une

| # | Difficulté | Origine mathématique |
|---|------------|----------------------|
| D1 | Alias entremêlés | chaque $\mu_k$ engendre un peigne $\mu_k+jm\omega$ ; $n_x$ peignes se superposent dans la zone de Brillouin |
| D2 | Parasites de troncature | les alias de bord $|m|\to h$ perdent en précision (le couplage vers $|n|>h$ est coupé) et dérivent en « queues de comète » de faible densité |
| D3 | Géométrie circulaire | l'axe imaginaire est replié modulo $\omega$ ; une distance euclidienne plate déchire un cluster centré près du bord $\pm\omega/2$ |
| D4 | Multiplicité coïncidente | des exposants répétés (systèmes découplés, symétriques) font coïncider exactement les peignes |
| D5 | Éclatement défectif | un exposant défectif de multiplicité $m$ se disperse en $\varepsilon^{1/m}$ (section 4.3) |
| D6 | Dynamique des multiplicateurs | $|\rho_k| = e^{T\operatorname{Re}\mu_k}$ couvre plusieurs décades, ce qui ruine tout regroupement dans l'espace des multiplicateurs |
| D7 | Cohérence de la période | la pulsation $\omega$ utilisée pour construire $H_h$ et pour replier doit coïncider ; un écart disperse tout le spectre |

Les difficultés D3, D4 et D6 ont été les pièges récurrents des tentatives antérieures, détaillés en section 7.

---

## 6. Algorithme retenu

L'algorithme `findTruelmV3` traite le spectre brut $\{\lambda_i\}$ ($N=(2h+1)n_x$ valeurs propres de $H_h$) en cinq étapes, toutes fondées sur une **métrique cylindrique unique**.

### 6.1 Repliement et métrique cylindrique

On replie chaque partie imaginaire dans la zone de Brillouin,

$$
\varphi_i = \operatorname{wrap}_\omega\!\big(\operatorname{Im}\lambda_i\big), \qquad \operatorname{wrap}_\omega(x) = \big((x+\tfrac{\omega}{2}) \bmod \omega\big) - \tfrac{\omega}{2},
$$

et on mesure les distances avec

$$
d^2(\lambda_a,\lambda_b) = \big(\operatorname{Re}\lambda_a-\operatorname{Re}\lambda_b\big)^2 + \operatorname{wrap}_\omega\!\big(\operatorname{Im}\lambda_a-\operatorname{Im}\lambda_b\big)^2 .
$$

Les exposants vivant dans $\mathbb{C}$, leurs parties réelle et imaginaire portent la même unité physique ($\mathrm{s}^{-1}$) : la métrique reste **brute, sans normalisation d'axe**. Diviser l'axe imaginaire par $\omega$ (comme le faisaient les versions antérieures) détruit cette géométrie et brise l'invariance d'échelle — c'est la cause directe de D3 mal traitée.

### 6.2 Densité cylindrique et seuil

On estime une densité à noyau gaussien en métrique cylindrique,

$$
\rho(\lambda_i) = \sum_{j=1}^{N} \exp\!\left(-\frac{d^2(\lambda_i,\lambda_j)}{2\,\mathrm{bw}^2}\right).
$$

Les alias propres d'un même mode coïncident après repliement et s'additionnent : le cœur d'un peigne atteint une densité élevée, tandis que les parasites de bord (D2), isolés, restent en bas. Un seuil $\tau$ écarte les pics parasites :

$$
\text{candidat pic} \iff \rho(\lambda_i) > \tau \cdot \max_j \rho(\lambda_j), \qquad \tau = 0.15 .
$$

### 6.3 Extraction gloutonne des pics distincts

On itère : sélection du point de densité maximale restant, puis suppression de son voisinage avec des **tolérances par axe**

$$
|\operatorname{Re}\lambda - \operatorname{Re}\lambda_{\text{pic}}| \le R_{\mathrm{re}}, \qquad \big|\operatorname{wrap}_\omega(\varphi-\varphi_{\text{pic}})\big| \le R_{\mathrm{im}} .
$$

Le rayon réel $R_{\mathrm{re}}$ est fixé par l'échelle réelle du spectre, le rayon imaginaire est plafonné par une fraction de la période, $R_{\mathrm{im}} = \min(R_{\mathrm{re}},\,0.15\,\omega)$, afin de ne pas sur-fusionner les zones de Brillouin étroites (petit $\omega$). La séparation des deux rayons traite D3 sans coupler la résolution réelle à $\omega$.

### 6.4 Affectation de Voronoï et multiplicité par la masse

Chaque point du spectre est affecté à son pic le plus proche (métrique cylindrique). La masse d'un pic est son nombre de points affectés, $\mathrm{mass}_k$. La multiplicité se lit alors par comparaison à la masse de référence d'un mode simple, $\rho_{\mathrm{ref}} = (2h+1) = N/n_x$ :

$$
m_k = \operatorname{round}\!\left(\frac{\mathrm{mass}_k}{\rho_{\mathrm{ref}}}\right).
$$

L'exposant $\mu_k$ (lu au point de densité maximale du cluster) est alors **répliqué $m_k$ fois**. C'est le point décisif pour D4 et D5 : un mode coïncident ou défectif de multiplicité $m$ ne se sépare pas — sa valeur est unique, seule sa multiplicité compte. Toute tentative de scission (par exemple un $k$-means local) d'un nuage de points coïncidents renvoie une copie correcte plus des artefacts.

### 6.5 Validation de Liouville

La formule de Liouville fournit une condition nécessaire exacte : la somme des exposants de Floquet égale la trace de l'harmonique continue $A_0$ de $A(t)$,

$$
\sum_{k=1}^{n_x} \mu_k = \operatorname{tr} A_0 .
$$

Le résidu relatif $\big|\sum_k\mu_k - \operatorname{tr}A_0\big| / |\operatorname{tr}A_0|$ sert de détecteur d'échec : un regroupement erroné le fait diverger. Sur le cas de référence vérifié, ce résidu vaut $2.4\times10^{-4}$.

**Limite à connaître :** la somme $\sum_k\mu_k$ porte sur les représentants repliés. Un mode situé exactement sur le bord $\pm\omega/2$ peut être représenté par $\mu$ ou $\mu\mp j\omega$, ce qui décale la somme d'un multiple de $\omega$ et gonfle artificiellement le résidu même si les exposants sont corrects. Le critère de Liouville n'est donc pas fiable sur le bord ; le pilote (section 8) utilise alors la stabilité en $h$.

---

## 7. Pourquoi les approches antérieures échouaient

### 7.1 Exponentiation vers les multiplicateurs (D6)

L'idée séduisante consiste à passer aux multiplicateurs $\rho_k = e^{\mu_k T}$, ce qui replie automatiquement chaque peigne sur un point unique, puisque $e^{(\mu_k+jm\omega)T}=\rho_k$. Aucun repliement manuel n'est nécessaire.

Le défaut est rédhibitoire. Le module $|\rho_k| = e^{T\operatorname{Re}\mu_k}$ explose en dynamique : sur le cas de référence ($T=2\pi$), un mode stable de partie réelle $-2{,}02$ donne $|\rho|\approx 3\times10^{-6}$, un mode instable de partie réelle $2{,}65$ donne $|\rho|\approx 1.6\times10^{7}$, soit un rapport de l'ordre de $10^{13}$. Tous les modes stables s'agglutinent indistinctement à l'origine, le conditionnement s'effondre, et l'arithmétique flottante perd les petits modules par soupassement. On a échangé un problème additif circulaire, réparable, contre un problème multiplicatif de dynamique, non réparable.

### 7.2 Clustering euclidien plat (D3)

Les solveurs `dbscan` et `kmeans` des versions antérieures regroupent sur la seule partie imaginaire repliée, en distance euclidienne plate. Deux défauts en résultent : ils ignorent la partie réelle (deux modes de même fréquence mais d'amortissement différent deviennent inséparables), et ils déchirent les clusters centrés près du bord $\pm\omega/2$ faute de métrique circulaire. Le solveur `robust` calculait bien une densité cylindrique, mais empilait ensuite un DBSCAN euclidien sur des coordonnées normalisées, réintroduisant le même biais de bord.

### 7.3 Scission d'un super-cluster coïncident (D4)

Le solveur `robust` détectait correctement la multiplicité par comptage ($D = \mathrm{mass}/\rho_{\mathrm{ref}}$), mais tentait ensuite de **scinder** le super-cluster par un $k$-means local en $m$ sous-modes. Sur un nuage de points rigoureusement coïncidents, aucune séparation n'existe : le $k$-means renvoie un centre correct et $m-1$ centres parasites. La réplication (section 6.4) remplace la scission et résout le cas.

---

## 8. Pilote de raffinement en $h$

Le coût d'une extraction est négligeable (de l'ordre de la milliseconde) devant celui du calcul du spectre $H_h$ (la décomposition propre, $O(N^3)$). La stratégie naturelle consiste donc à augmenter la troncature jusqu'à convergence :

$$
h \leftarrow \lceil 1.5\,h \rceil, \qquad \text{tant que } \Delta\mu(h) \ge \mathrm{tol},
$$

où $\Delta\mu(h)$ est l'écart cylindrique maximal entre les exposants obtenus à deux valeurs successives de $h$, après appariement unique. Le critère de stabilité en $h$ est préféré au résidu de Liouville parce qu'il reste valide sur le bord de Brillouin (section 6.5). La convergence exige deux $h$ consécutifs concordants, ce qui empêche un $h$ fortuit de déclencher un arrêt prématuré.

Ce mécanisme matérialise la division des rôles : `findTruelm` n'a pas à résoudre les faibles troncatures, c'est le pilote qui redemande un meilleur spectre à $H_h$. **La seule exigence est la robustesse aux grands $h$, toujours, y compris en présence de modes dégénérés.**

---

## 9. Validation

### 9.1 Cas régulier et multiplicités

Sur le cas de référence (matrice $3\times3$ périodique, $h$ stocké $=15$) et ses constructions dégénérées, l'extraction est exacte pour $h\ge 10$ :

| Cas | Multiplicités | Résultat |
|-----|---------------|----------|
| Référence (paire conjuguée partageant $\operatorname{Re}$) | $1,1,1$ | exact |
| `blkdiag(A,A)` | $2,2,2$ | exact |
| `blkdiag(A,A,A)` | $3,3,3$ | exact |
| Mixte $A,A,$ scalaire | $2,2,2,1$ | exact |

La robustesse se maintient jusqu'à $h=160$ sans écart (le coût restant dominé par la décomposition propre, l'extraction demeurant à $\sim 60\ \mathrm{ms}$ pour $N=2889$).

### 9.2 Cas défectifs

Construits par conjugaison d'un bloc de Jordan, $A = S\,J_m(-1)\,S^{-1}$, puis rendus **périodiques** par la construction de Floquet-Lyapunov $A(t) = \dot{R}(t)R(t)^{-1} + R(t)\,J\,R(t)^{-1}$ avec $R$ une rotation périodique. Cette construction garantit des exposants exactement égaux à $\operatorname{eig}(J)=-1$ (défectif), tout en produisant un système périodique plein dont la troncature génère une réelle pollution. L'extraction renvoie $-1$ avec la bonne multiplicité, le résidu suivant l'éclatement $\varepsilon^{1/m}$ : de $3\times10^{-9}$ ($m=2$) à $4\times10^{-2}$ ($m=10$).

Un constat important guide le choix des cas de test : une matrice **constante** produit une matrice de Hill exactement bloc-diagonale, donc des peignes propres sans couplage inter-blocs ni pollution de troncature. Les cas dégénérés constants sous-estiment la difficulté ; seuls les systèmes réellement périodiques exercent la robustesse visée.

### 9.3 Plancher de résolution

Deux modes distincts plus proches que la tolérance de fusion sont irrémédiablement confondus à une troncature donnée. Le balayage synthétique le situe vers un écart de $0.04$ (en unités de l'axe, métrique cylindrique). Ce comportement est **détecté**, non masqué : le résidu de Liouville ou la non-convergence du pilote signalent l'échec, qui se corrige en augmentant $h$ (les peignes se resserrent). Aucun échec n'est silencieux.

---

## 10. Limites ouvertes

- **Coalescence sous le plancher de résolution.** Deux exposants séparés de moins de $\sim 0.04$ fusionnent ; seul un $h$ plus grand, ou une information supplémentaire, les sépare.
- **Modes sur le bord $\pm\omega/2$.** Le résidu de Liouville y est inexploitable (ambiguïté de branche $m\omega$) ; le critère de stabilité en $h$ le remplace, mais une caractérisation propre du bord reste à formaliser.
- **Blocs de Jordan d'ordre élevé.** L'éclatement $\varepsilon^{1/m}$ croît avec $m$ ; au-delà de $m\approx 8$–$10$ sur un Hill mal conditionné, l'éclatement peut dépasser la tolérance de fusion et fragmenter le comptage de multiplicité.
- **Vraie dégénérescence de représentation.** Deux exposants distincts de même partie réelle dont les parties imaginaires diffèrent d'un multiple exact de $\omega$ coïncident sur le cylindre : ils sont alors le même exposant de Floquet (multiplicité $2$), et la distinction n'a pas de sens sur le seul spectre. Une séparation exigerait l'information des vecteurs propres du problème de Hill.

---

## 11. Fichiers

| Fichier | Rôle |
|---------|------|
| `findTruelmV3.m` | solveur d'extraction (sections 6.1–6.5) |
| `findTrueFloquetV3.m` | pilote de raffinement en $h$ (section 8) |
| `test_findTruelmV3.m` | tests de régression sur spectres `HmqNEig` réels |
| `stress_v3.m` | suite synthétique des cas limites (sections 5, 9.3) |
| `fixture_floquet_hardcase.mat` | cas de référence archivé (matrice, exposants vrais, $\operatorname{tr}A_0$) |
