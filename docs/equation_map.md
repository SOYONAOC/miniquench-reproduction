# 方程、单位、实现与验证

原文固定为 arXiv:2609.19265v1。以下首列为 HTML 编号；PDF 第 2 节对应 2.1–2.21，第 3 节另行编号，不能将 HTML 编号直接作为 PDF 编号引用。全文与图注已核对 PDF、HTML、TeX，原始文件和 SHA256 见 data/manifest.json。

|原文|定义与单位|实现|验证|
|---|---|---|---|
|§1 宇宙学|Ωm=.3111, ΩΛ=.6889, Ωb=.0489, h=.6766；平直、H(z) s⁻¹|physics.hubble|不调用 AuroraLF 或 Astropy 预设宇宙学|
|(1),(2),(21),(22)|Menc(r)=2σ²r/G, σ²=GMh/rvir；Mh 是总质量尺度，Menc 是等温球包围质量|physics.enclosed_mass、dynamics.shell_budget|引力随 r⁻¹；factor=1 仅为标记分支|
|(3)|rvir=[2GMh/(Δh Ωm H²)]^(1/3)，Δh=200，cm|physics.halo|逐项保留 Ωm 与 H(z)；不换通用 virial 定义|
|(4)|dMh/dt=A Mh (1+z)^(5/2)，A=.03/Gyr|physics.growth_rate、halo|解析指数增长；Myr 内部时钟|
|(5),(6)|dMg/dt=fb dMh/dt−SFR；SFR=fstar Mg/tff，tff=.0026/H|physics.reservoir|解析气体解和 Mg+Mformed 质量恒等式|
|§2.2|fstar 为每自由落体效率；ftilde 只映射观测总恒星质量|reproduce.plot_results|本轮 formed_msun 独立保存；映射取 ftilde=.01,.1|
|§2.3|R0=κ rvir(t)，κ=.035|physics.boundary|移动边界速度、加速度，事件用 r−R0|
|(7)|Φ=dN/dm，∫mΦdm=1，每形成 1 Msun；m=.08–100，β=1.6|feedback.IMF|五组斜率/特征质量的质量积分|
|(8)|τ=10⁴ m^(−2.5) Myr；SN 前身星 8–100 Msun|feedback.lifetime、sn_kernel|核积分=∫8^100 Φdm；核支持 .1–55.2427 Myr|
|(9),(10)|UV 黑体波段 100–4000 Å，Teff=5772 m^(5/8) K；LUV=∫SFH lUV dage|feedback.Kernels、History|瞬时、恒定和关停后响应；独立积分比较|
|§2.6.1 半径缺项|R*=Rsun m^0.5 为明确推断，非已确认作者公式|Config.stellar_radius_power|0.7 分支；不能据此宣称 exact reproduction|
|(11)|FUV=(1−exp(−τUV)+τIR)LUV/c；κUV=10³, κIR=10^.7 cm²/g，D=Zsun_fraction/162|radiation.radiation_forces|CGS 力 dyne；质量为覆盖部分，面积同乘 fcover|
|(12)–(14)|lLyα=.9×(2/3) Qion hc/1215.67Å；FLyα=MF LLyα/c|feedback.Kernels、radiation.radiation_forces|关停后仍有光子；Nebrin 原始 MF 函数|
|Nebrin et al. 2025 §2.3|含尘埃吸收、反冲、光子破坏、可选速度梯度的 MF 拟合|third_party/Lyman-alpha-feedback/M_F_fit|固定 commit 与校验值；不以常数代替|
|Schure et al. 2009 表2、式3|Lcool=Λhd(T)nH²V，Λhd 已含 ne/nH|data/schure_cooling.csv、radiation.thermal_equilibrium|110 行精确解析；不重复乘 ne/nH；越界/无稳定根显式报错|
|(15)|SN rate=∫SFR(t′)Φ(m(age)) |dm/dage|dt′，SN/Myr|feedback.History|核有限范围、top-hat 独立积分、因果性|
|(16)–(18)|Pdot=LSN/(2πr³)−P/tcomp−5Pv/r；Fp=4πr²P；LSN=εk Ndot 10⁵¹ erg|physics.pressure_rhs、dynamics|Pr⁵ 守恒；固定 r 指数冷却|
|(17)|tcomp=120[(1+z)/10]⁻⁴ Myr|physics.pressure_rhs|冷却时间一个 e-fold 精度测试|
|(19),(20)|fcover=Nfil θfil/(4π)，Nfil=3；扫描直接指定覆盖率|Config.cover|fcover=1 恢复全壳层；吸积率不随 fcover 下降|
|(21)|mseg vdot=−G Menc mseg/r²−mdotseg(vin+v)+fcover ΣFout|dynamics.shell_budget、run|d(mv)/dt=ΣFexternal−mdot vin；回落分支；总气体守恒|
|§2,图2|tSF=离开形成区；tquench=返回−启动|dynamics.run|转向不结束淬火；不返回只有下限；窗口结束不强制 SFR=0|
|(25)–(33)|解析 Mmax=4.4e11(fstar/.01)^1.5 × 分段时间函数|analytic.analytic_mass|z=6；55 Myr 接合连续；独立于动态轨道|
|§3.3 图4|仅 SN+引力，Nsn 数值卷积和积分，Fp=2Nsn E/R0|analytic.numerical_mass|随质量缩放解力平衡；增长和静态分开；无辐射/冷却/动量项|
|(34),(35)|Mstar,total≈ftilde fb Mh|reproduce.plot_results|不得拿本轮形成质量替换观测总质量|

## 积分约定

内部时间 Myr，存储长度 cm、速度 km/s、压力 erg/cm³、质量 Msun、SFR Msun/yr、力 dyne。求导使用无量纲半径/压力/质量状态并显式转换 Myr↔s。历史格距、RK45 求解器最大步长、输出间隔分别配置。气体储库的解析形成历史是两个指数项之和；历史卷积预积分 exp(−rate×age) 加权年龄核，在 [max(0,t−tSF),t] 取积分差。指数权重在年龄格中心近似，核积分来自累积 IMF 表；通过 history_dt 收敛检验。非网格关停时刻的反馈连续。RHS 不写入历史。

## 资料边界

全模型热学闭合未由论文唯一确定，见 ambiguities.md。`full_approx` 含所有反馈项但不是作者实现的逐行复原；`sn_only` 是明确的调试对照。无作者数值数据，参考曲线来自其源包矢量 PDF 的路径坐标，具有绘图降采样和读图误差。
