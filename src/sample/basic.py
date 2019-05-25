#!/usr/bin/env python3

# This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

import math
from consts import Consts


class Cell():
    def __init__(self, id=None, pos=[0, 0], veloc=[0, 0], radius=5):
        # ID to judge Player or free particle
        self.id = id
        # Variables to hold current position
        self.pos = pos
        # Variables to hold current velocity
        self.veloc = veloc
        # Variables to hold size
        self.radius = radius

        # Properties
        self.collide_group = None
        self.dead = False

    # Methods
    def distance_from(self, other):
        """Calculate the distance from another cell.

        Args:
            other: another cell.
        Returns:
            the minimum distance.

        """
        dx = self.pos[0] - other.pos[0]
        dy = self.pos[1] - other.pos[1]
        min_x = min(abs(dx), abs(dx + Consts["WORLD_X"]), abs(dx - Consts["WORLD_X"]))
        min_y = min(abs(dy), abs(dy + Consts["WORLD_Y"]), abs(dy - Consts["WORLD_Y"]))
        return (min_x ** 2 + min_y ** 2) ** 0.5

    def collide(self, other):
        """Determine if it collides with another cell.

        Args:
            other: another cell.
        Returns:
            True / False.

        """
        return self.distance_from(other) < self.radius + other.radius

    def area(self):
        """Calculate the area of the cell.

        Args:

        Returns:
            the area of the cell.

        """
        return math.pi * self.radius * self.radius

    def stay_in_bounds(self):
        """Make the out-of-bounds cell stay within the bounds.

        Args:

        Returns:


        """
        if self.pos[0] < 0:
            self.pos[0] += Consts["WORLD_X"]
        elif self.pos[0] > Consts["WORLD_X"]:
            self.pos[0] -= Consts["WORLD_X"]

        if self.pos[1] < 0:
            self.pos[1] += Consts["WORLD_Y"]
        elif self.pos[1] > Consts["WORLD_Y"]:
            self.pos[1] -= Consts["WORLD_Y"]

    def limit_speed(self):
        """Enforce speed limits.

        Args:

        Returns:


        """
        if self.veloc[0] > Consts["MAX_VELOC"]:
            self.veloc[0] = Consts["MAX_VELOC"]
        elif self.veloc[0] < -Consts["MAX_VELOC"]:
            self.veloc[0] = -Consts["MAX_VELOC"]

        if self.veloc[1] > Consts["MAX_VELOC"]:
            self.veloc[1] = Consts["MAX_VELOC"]
        elif self.veloc[1] < -Consts["MAX_VELOC"]:
            self.veloc[1] = -Consts["MAX_VELOC"]

    def move(self, frame_delta):
        """Move the cell according to its velocity.

        Args:
            frame_delta: Time interval between two frames.
        Returns:


        """
        self.collide_group = None
        # Adjust the position, according to velocity.
        self.pos[0] += self.veloc[0] * frame_delta
        self.pos[1] += self.veloc[1] * frame_delta
        self.stay_in_bounds()
        self.limit_speed()


# 地图是1000*500的，但是由于边界的穿越性，本来的反弹变成了直接到对面
# 那就相当于将 地图扩展成了3*3个1000*500
# 计算距离和其他参数的时候，我们只需要考虑最近的那一种情况就行了

# -------------------------------------------
class Player():
    search_radius = 200  # (搜索半径）
    danger_radius = 50  # 危险判定半径
    vmax = 0.7  # 最高的速度，再快就不追了，自动捕获即可

    def __init__(self, id, arg=None):
        self.id = id

    # -------------------------------辅助计算函数----------------------------------

    # 9个（包括本体）镜像的球
    def mirro_cell(self, cell):
        mirror = []
        x = cell.pos[0]
        y = cell.pos[1]
        for i in range(-1, 2):
            for j in range(-1, 2):
                temp = Cell(cell.id, [0, 0], cell.veloc, cell.radius)
                temp.pos[0] = x + i * Consts["WORLD_X"]
                temp.pos[1] = y + j * Consts["WORLD_Y"]
                mirror.append(temp)
        return mirror

    # b相对a的径向速度（仅考虑未扩展的地图）
    def delta_vr(self, a, b):
        dx = -a.pos[0] + b.pos[0]
        dy = -a.pos[1] + b.pos[1]
        dvx = -a.veloc[0] + b.veloc[0]
        dvy = -a.veloc[1] + b.veloc[1]
        dv = (dvx ** 2 + dvy ** 2) ** 0.5
        if dvx == dvy == 0:
            return 0
        else:
            dvr = dvx * (dx / dv) + dvy * (dy / dv)
            return dvr

    # b相对a的速度大小
    def delta_v(self, a, b):
        dvx = -a.veloc[0] + b.veloc[0]
        dvy = -a.veloc[1] + b.veloc[1]
        dv = (dvx ** 2 + dvy ** 2) ** 0.5
        return dv

    # 直接计算两个点之间的距离
    def delta_r(self, a, b):
        dx = -a.pos[0] + b.pos[0]
        dy = -a.pos[1] + b.pos[1]
        return (dx ** 2 + dy ** 2) ** 0.5

    # 匀速运动的a，b之间的最小距离(a到b所在直线的最小距离）
    def min_distance(self, a, b):
        dx = -a.pos[0] + b.pos[0]
        dy = -a.pos[1] + b.pos[1]
        dvx = -a.veloc[0] + b.veloc[0]
        dvy = -a.veloc[1] + b.veloc[1]
        dr = (dx ** 2 + dy ** 2) ** 0.5
        theta1 = math.atan2(dvy, dvx)  # 相对速度与x轴夹角
        theta2 = math.atan2(dy, dx)  # 相对位置于x轴 夹角
        return abs(dr * math.sin(theta2 - theta1))

    # 匀速运动到最小距离所花费的时间
    def min_time(self, a, b):
        return self.min_distance(a, b) / self.delta_v(a, b)

    def sumR(self, a, b):
        return a.radius + b.radius

    # ------------------------------策略函数----------------------------------

    # a到b的最短迎接方案（b向着a过来）
    def best_method(self, a, b):
        """
        :type a: Cell
        :type b: Cell
        """
        # 先不考虑中间有球
        mirror = self.mirro_cell(b)
        ways = []
        if a.distance_from(b) - self.sumR(a, b) > self.search_radius:
            return None
        #  分析漂浮到达每一个镜像的情况，计算到达二者相距最小距离的时间和这个最小的方向
        for i in mirror:
            if self.delta_r(a, i) - self.sumR(a, i) > self.search_radius:
                continue
            dvr = self.delta_vr(a, i)
            if dvr >= 0:  # 将反向的排除暂时
                continue
            min_dis = self.min_distance(a, i)
            R = self.sumR(a, i)
            dr = self.delta_r(a, i)
            dv = self.delta_v(a, i)
            if min_dis <= R:  # 这个是可以自己碰撞到的，
                time = (math.sqrt(dr ** 2 - min_dis ** 2) - math.sqrt(R ** 2 - min_dis ** 2)) / dv
                if a.veloc[0] ** 2 + a.veloc[1] ** 2 >= self.vmax ** 2:  # 上限速度
                    ways.append((time, None))
                else:
                    dx = -a.pos[0] + i.pos[0]
                    dy = -a.pos[1] + i.pos[1]
                    theta = math.atan2(dx, dy) + math.pi  # 极轴是y轴,反向喷射球,哪个方向加速可以帮它最快到达？？目前假设径向
                    ways.append((time, theta))
            else:  # 主动凑上去
                time = min_dis / dv
                if a.veloc[0] ** 2 + a.veloc[1] ** 2 >= self.vmax ** 2:
                    ways.append((time, None))
                else:
                    dvx = -a.veloc[0] + i.veloc[0]
                    dvy = -a.veloc[1] + i.veloc[1]
                    theta = -math.atan2(dvy, dvx) + math.pi
                    ways.append((time, theta))
        if not ways:
            return None
        best = sorted(ways, key=lambda way: way[0])[0]  # 最短时间的选择！
        return best

    # a如何逃离直接背b吃掉的命运
    def best_escape(self, a, b):
        dangers = []
        if a.distance_from(b) - self.sumR(a, b) > self.danger_radius:
            return None
        mirror = self.mirro_cell(b)
        for i in mirror:
            dvr = self.delta_vr(a, i)
            if dvr >= 0:  # 将反向的排除暂时
                continue
            R = self.sumR(a, i)
            dr = self.delta_r(a, i)
            if dr - R > self.danger_radius:
                continue
            min_dis = self.min_distance(a, i)
            dv = self.delta_v(a, i)
            if min_dis <= R:  # 这个是可以自己碰撞到的,有危险 (在加大一点？不要极限操作吧...）
                time = (math.sqrt(dr ** 2 - min_dis ** 2) - math.sqrt(R ** 2 - min_dis ** 2)) / dv
                dx = -a.pos[0] + i.pos[0]
                dy = -a.pos[1] + i.pos[1]
                theta = math.atan2(dx, dy)  # 极轴是y轴,喷射球,指向最近距离方向
                dangers.append((time, theta))  # 返回碰撞时间和逃离角度
        if not dangers:
            return None
        most_danger = sorted(dangers, key=lambda danger: danger[0])[0]
        return most_danger

    # 在给出的会造成危险的情况下，多个危险怎么跑？，限定到只有两个？
    def escape(self, dangers):
        if len(dangers) == 1:
            print('danger!')
            return dangers[0][1]
        if len(dangers) >= 2:
            dangers = sorted(dangers, key=lambda danger: danger[0])[:2]
            a = dangers[0]
            b = dangers[1]
            theta = (a[1] * b[0] + b[1] * a[0]) / (a[0] + b[0])  # 角度加权紧急程度
            print('danger!')
            return theta
        if len(dangers) == 0:
            return None

    def strategy(self, allcells):
        player: Cell = allcells[self.id]
        enemy = None
        bigger = []
        smaller = []
        # 标记出敌人，区分大的和小的
        for i in allcells:
            if i == player:
                continue
            if i.id + self.id == 1:  # 只能是0+1了
                enemy = i
            if i.radius >= player.radius:
                bigger.append(i)
            else:
                smaller.append(i)

        # 逃跑：
        if bigger == []:
            pass
        else:
            dangers = []
            for cell in bigger:
                danger = self.best_escape(player, cell)
                if danger is not None:
                    dangers.append(danger)
            if not dangers:
                pass
            else:
                return self.escape(dangers)

        # 吃球
        if len(smaller) == 0:
            return None
        methods = []
        for cell in smaller:
            method = self.best_method(player, cell)
            if method is not None:
                methods.append(method)
        if not methods:
            return None
        best_method = sorted(methods, key=lambda m: m[0])[0]
        # if best_target.radius < player.radius * 0.11:  # 防止去吃太小的
        #     return None

        v = (player.veloc[0] ** 2 + player.veloc[1] ** 2) ** 0.5
        if v > self.vmax:
            return None
        else:
            print('eat!')
            return best_method[1]
