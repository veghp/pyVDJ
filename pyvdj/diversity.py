# Copyright 2019 Peter Vegh
#
# This file is part of pyVDJ.
#
# pyVDJ is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyVDJ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with pyVDJ.  If not, see <https://www.gnu.org/licenses/>.


from math import log


def shannon(data):
    # Shannon diversity index
    # 'data': dictionary of species: count
    N = sum(data.values())
    summa = [ (float(n)/N) * log(float(n)/N) for n in data.values() if n != 0 ]

    return -sum(summa)


def simpson(data):
    # Simpson diversity index
    # 'data': dictionary of species: count
    N = sum(data.values())
    summa = [ (float(n)/N)**2 for n in data.values() if n != 0 ]

    return sum(summa)
