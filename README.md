# asti 
This is cagd library with the following functionality implemented so far
- splines (regular,periodic,rational).
- conic arc, circle, line.
- subdivide curve.
- split splines to bezier segments.
- extend splines.
- interpolate points to form splines (both periodic and regular) using pchip.
- raise degree.
- insert knots.
- remove knot for curve fairing.
- cubic and quad approximation
- distance of a point to a curve.
- arc lengths of curves.
- implicitize 2d curves
- monomial and Bezier-Bernstein form of polynomial curves and change of basis.
- point inversion on spline to determine parameters.
- bspline surface, ruled surface, extrude.
[![Build Status](https://travis-ci.org/svark/asti.svg?branch=master)](https://travis-ci.org/svark/asti)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/svark/asti?branch=master&svg=true)](https://ci.appveyor.com/project/svark/asti)
<br/>
<a href="https://scan.coverity.com/projects/5900">
  <img alt="Coverity Scan Build Status"
       src="https://scan.coverity.com/projects/5900/badge.svg"/>
</a>

My plan is to read up text books Hoschek/Farin/Piegl and translate the book ideas into code as much as I can in here.
The project uses biicode to handle builds and catch to perform unit testing.



