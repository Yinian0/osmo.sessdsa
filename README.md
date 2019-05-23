# osmo.sessdsa

地空数算课2019期末作业：星际吞噬

## 游戏概述

本游戏的原型为[Osmos](https://www.osmos-game.com/)，一个由Hemisphere Games开发的独立游戏，曾在2009年[独立游戏节](https://zh.wikipedia.org/独立游戏节)获得3个奖项的提名。原作支持多个平台，欢迎通过官方渠道购买/下载。

Osmo是受到到Osmos启发，而制作的一个简化版游戏框架。这里有一个部署在Heroku上的JavaScript单机版本，可以[在线体验](https://osmoss.herokuapp.com/)。

## 规则介绍

本次大作业的内容是编写一个AI函数，来进行Osmo游戏，并取得胜利。具体来说，大作业采用的是双玩家对战。AI要做的事情，就是使之作为『玩家』做出决策控制自身星体，不断成长并通过『**吞噬**』（吸收）由另一个AI操控的星体获胜。

为了保证比赛公平性，并适当控制AI编写难度（亦可以使人类玩家获得更好的游戏体验），下面将详细介绍大作业所采用的规则。其它的一些设计思路可以参见[IDEAS.md](IDEAS.md)。

### 游戏基本内容

游戏场地限制为`1000x500`的矩形。游戏内基本单元为二维平面上的一些球（或者说，圆盘），称为『**星体**』。每个星体有一定的**半径**，半径决定了其面积，并定义星体的『**质量**』正比于面积。星体互不重叠且可以一定**速度**运动。

### 星体编号与玩家

游戏具有两个玩家，每个玩家将控制一个星体。如果由玩家控制的星体被吞噬（见『[接触与吞噬判定](#接触与吞噬判定)』部分），则这名玩家失利。

对局开始时，随机确定一名『**主玩家**』，使其控制编号0的星体，另一位玩家控制编号为1的星体。两玩家星体分别位于坐标`(250, 250)`与`(750, 250)`处，初始半径均为`DEFAULT_RADIUS`，初始速度静止。  
场地中除了玩家星体外还将随机分布`CELLS_COUNT`个星体，其半径分为大、中、小三类，初始速度范围为-1至1。随机分布的星体将获得递增的编号。在星体被吞噬后，标记其死亡并不再更新；每次在弹射产生新的星体后，使其获得一个新的编号。

换句话说，两玩家星体编号分别为0、1，编号2及以后为非玩家星体。

### 漂浮与边界

在没有外界影响的情况下，所有星体将做匀速直线运动；而当星体部分穿出游戏场地的边界时，穿出的部分将速度不变地移动至场地相对的一侧。

### 接触与吞噬判定

当两个或多个星体相接触（即二者中心距离小于二者半径之和）时，将会产生『**吞噬**』，吞噬过程保持『质量』守恒和『**动量**』守恒。被吞噬的星体将会从对局中移除。  

对于二体接触，按以下规则判定何者会吞噬/被吞噬：
- 当二者半径不等时，半径较大者吞噬半径较小者。
- 当二者半径在误差范围（根据Python使用的IEEE-754浮点数精度）内相等，星体编号较小者吞噬星体编号较大者。

当一个星体同时与多个星体接触，或多个星体同时相互接触（『**同时**』指一个运算帧内）时，按以下规则处理：
- 对于相接触的一组星体，由其中半径最大者吞噬其它全部星体。
- 半径最大者多于一个时，由其中星体编号最小者吞噬其它全部星体。

若一次吞噬发生后又产生了新的接触情况，则新的接触情况在下一帧处理（可能下一帧已脱离接触，则不再处理）。

### 弹射

玩家控制的星体除遵守以上规则以外，还可主动将自身的部分质量作为一个新的星体『**弹射**』出去。弹射角度可由玩家自由控制。弹射过程视为碰撞的逆过程，依然遵守『质量』和『动量』守恒。  
玩家可以利用在弹射星体过程中的反作用来控制自身运动，进而不断吞噬比自身小的星体，并躲避比自己更大的星体。

玩家每一帧中最多可弹射1次星体，弹射星体的质量与本体质量之比为一恒定值`EJECT_MASS_RATIO`，相对自身速度为恒定值`DELTA_VELOC`。

## 游戏结束判定

后文中的『代码出错』定义为：抛出未被捕获的异常，或者返回值不满足输入输出要求。

### 获胜条件

- 对方星体被吞噬且己方存活
- 达到最大回合数且己方星体半径大于对方
- 对方代码出错且己方运行正常

### 失利条件

- 己方星体被吞噬且对方存活
- 达到最大回合数且对方星体半径大于己方
- 己方代码出错且对方运行正常

### 平局条件

- 双方星体同时被吞噬
- 达到最大回合数且双方星体半径相等
- 双方代码同时出错

## 数据结构与算法

对于游戏内核的描述可以参见[KERNEL.md](KERNEL.md)。

## 时间空间限制

双方代码的运算总时长均限制为10s。在运算时长耗尽后，玩家将会失去对星体的控制（只是无法再进行弹射，而非立即判负），直到游戏以某种方式结束。

## 备注

欢迎通过Issue或者Pull Request帮助我们完善代码、修复问题。

## License

Licensed under the GPLv3 license.
