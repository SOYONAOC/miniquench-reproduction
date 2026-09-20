# 逐图参数与采样

全局：固定 z；Ωm=.3111、ΩΛ=.6889、Ωb=.0489、h=.6766；κ=.035、εk=1、Z/Zsun=1e−4；初始 Mg=fb Mh（特别说明除外）；无先前恒星反馈、P=0、r=R0；IMF .08–100 Msun，β=1.6、正文 α=2.35；启动窗口 100 Myr。

|图|原图参数/范围|本项目实际采样|
|---|---|---|
|1|α=1.35/2.35；mchar=.2,5,10；.08–100 Msun|对数 1000 点；另画图5图注字面 α=−2.35；未伪造 Chabrier 原始经验函数|
|2|Mh0=10⁸,10⁹ Msun；fstar=.001；Chabrier-like；fcover=1|z=6 为上下文默认（图注未再次给出）；轨道到返回或 700 Myr；输出全部力/反馈/质量与事件|
|4|z=3,6,9；fstar=.01；mchar=.2,5；tSF=10–100 Myr|181 点；SN+引力简化；增长/固定暗晕对照；式33仅在 z6 叠加|
|5|z=6；fstar=.001,.01,.1；另固定 .01 扫 mchar=.2,1,5；Mh=10⁸–10¹³ Msun|quick 9 质量点，full 31 点；启动状态边界额外二分 4 次；ftilde=.01–.1 仅作横坐标映射|
|8|z=6；fstar=.001,.01；fcover=1,.1,.05,.005；Mh=10⁸–10¹³|同上；额外 fcover=1、初始 Mg=.1 fb Mh；吸积率保持原值|
|3|z=6，两质量，两 IMF 的反馈力比例|核心的作者曲线符合性未通过之前不宣称本图复现|
|6|形成时间/SFR 的解析质量边界|核心优先；完成状态见报告|
|7|z=3,6,9、fstar=.001、Chabrier-like|核心优先；100 Myr 为启动窗口|
|9|mini-quench 比例及 HMF|核心优先，未引入 Colossus|
|10|Compton 比例的形成阶段诊断|冷却单元测试已验证；逐图状态见报告|

## 工程参数（不是论文输入）

quick：history_dt=.2 Myr，max_step=1 Myr，output_dt=1 Myr。
full：history_dt=.05 Myr，max_step=.5 Myr，output_dt=.5 Myr。
rtol=2e−7，atol=1e−9（无量纲状态）；最大演化 700 Myr；最多 4 进程且每进程 BLAS 单线程。

这些配置及每一条实际运行配置均写入 configs/*.json 和 data/curves/summary_*.json。full_approx 的热学/恒星半径选择及有限分支详见 ambiguities.md；SN-only 不能替代完整反馈复现。
