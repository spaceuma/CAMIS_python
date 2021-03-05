"""
MIT License
-----------

Copyright (c) 2021 University of Malaga
Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use,
copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following
conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

Authors: J. Ricardo Sánchez Ibáñez, Carlos J. Pérez del Pulgar
Affiliation: Department of Systems Engineering and Automation
Space Robotics Lab (www.uma.es/space-robotics)
"""
# -*- coding: utf-8 -*-
"""
Simple test on inclined surface
"""

import numpy as np
from context import camis

x = np.linspace(0,50,51)
y = np.linspace(0,50,51)

XX,YY = np.meshgrid(x,y)

DEM = YY/5
demRes = 1.0
planRes = 1.0


# We create a new environment
env1 = camis.AnisotropicMap(DEM, demRes, planRes, (0,0))
env1.smoothMap(1.0)

cdRoots =  [0.0, 0.10, 1.0]
caRoots =  [0.0, 0.30, 1.0]
cl1Roots = [0.0, 0.20, 1.0]
cl2Roots = [0.0, 0.20, 1.0]

r1 = camis.CamisModel.fromRoots(cdRoots,caRoots,cl1Roots,cl2Roots, 45.0)
env1.computeVecCostMap(r1)

start = np.asarray([10,10])
goal = np.asarray([40,40])

env1.executeBiPlanning(goal,start)

#env1.executePlanning(goal,start)
#env1.showResults()
#
#fig, axes = plt.subplots(constrained_layout=True)
#env1.showMap('elevation',fig,axes)