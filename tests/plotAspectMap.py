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

plt.style.use('default')
plt.rcParams["font.family"] = "Constantia"
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
fig, ax1 = plt.subplots(figsize=(3.2, 3.9),nrows = 1, ncols = 1, constrained_layout=True)
cc = ax1.scatter(env.hexXmap, env.hexYmap, c = 180/np.pi*np.arctan2(env.hexAspectMap[1],env.hexAspectMap[0]),
                 cmap="gist_rainbow",s=16.0, vmin = -180.0, vmax = 180.0, rasterized=True)
cbar = fig.colorbar(cc, ax=ax1,shrink=0.6, location = 'top', ticks=[-180, -90, 0, 90, 180])
cbar.set_label('Aspect Angle [deg]')
ax1.set_xlim([0,48.0])
ax1.set_ylim([0,48.0])
ax1.set_xlabel('X [m]')
ax1.set_ylabel('Y [m]')
ax1.set_aspect('equal')
plt.savefig('numerical_tests_aspectMap_reduced.pdf',dpi=300)