r"""
Utilities for calculus and graphical outputs

This module defines helper functions which are used for simplifications of
symbolic expressions or for graphical outputs. 

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013-2015) : initial version

"""

#******************************************************************************
#       Copyright (C) 2015 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2015 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.sage_object import SageObject
from sage.version import version

def simple_determinant(aa):
    r"""
    Compute the determinant of a square matrix.

    This function, which is based on Laplace's cofactor expansion, is a
    workaround for a bug in Sage method
    :meth:`~sage.matrix.matrix2.Matrix.determinant`
    (cf. http://trac.sagemath.org/ticket/14403).

    NB: this bug was fixed in Sage 6.2.

    EXAMPLE::

        sage: a = matrix([[sqrt(x),0,0,0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        sage: a.determinant() # this resulted in an error in Sage < 6.2
        sqrt(x)
        sage: from sage.geometry.manifolds.utilities import simple_determinant
        sage: simple_determinant(a)
        sqrt(x)
    
    """
    from sage.matrix.constructor import matrix
    n = aa.nrows()
    if n == 1:
        return aa[0,0]
    res = 0
    sign = True
    for i in range(n):
        b = []
        for k in range(i):
            r = []
            for l in range(1,n):
               r.append(aa[k,l])
            b.append(r)
        for k in range(i+1,n):
            r = []
            for l in range(1,n):
               r.append(aa[k,l])
            b.append(r)
        bb = matrix(b)
        if sign:
            res += aa[i,0] * simple_determinant(bb)
        else:
            res -= aa[i,0] * simple_determinant(bb)
        sign = not sign
    return res

def simplify_sqrt_real(expr):
    r"""
    Simplify sqrt in symbolic expressions in the real domain.

    EXAMPLES:

    Simplifications of basic expressions::

        sage: from sage.geometry.manifolds.utilities import simplify_sqrt_real
        sage: assume(x<0)
        sage: simplify_sqrt_real( sqrt(x^2) )
        -x
        sage: simplify_sqrt_real( sqrt(x^2-2*x+1) )
        -x + 1
        sage: simplify_sqrt_real( sqrt(x^2) + sqrt(x^2-2*x+1) )
        -2*x + 1

    This improves over Sage's
    :meth:`~sage.symbolic.expression.Expression.canonicalize_radical` which yields
    incorrect results when x<0::

        sage: sqrt(x^2).canonicalize_radical() # wrong output
        x
        sage: sqrt(x^2-2*x+1).canonicalize_radical() # wrong output
        x - 1
        sage: ( sqrt(x^2) + sqrt(x^2-2*x+1) ).canonicalize_radical() # wrong output
        2*x - 1

    """
    from sage.symbolic.ring import SR
    from sage.calculus.calculus import maxima
    from sage.functions.other import sqrt
    # 1/ Search for the sqrt's in expr
    sexpr = str(expr)
    if 'sqrt(' not in sexpr:  # no sqrt to simplify
        return expr
    if 'D[' in sexpr:
        return expr    #!# the code below is not capable of simplifying
                       # expressions with symbolic derivatives denoted by Pynac
                       # symbols of the type D[0]
    pos_sqrts = []   # positions of the sqrt's in sexpr
    the_sqrts = []   # the sqrt sub-expressions in sexpr
    for pos in range(len(sexpr)):
        if sexpr[pos:pos+5] == 'sqrt(':
            pos_sqrts.append(pos)
            parenth = 1
            scan = pos+5
            while parenth != 0:
                if sexpr[scan] == '(': parenth += 1
                if sexpr[scan] == ')': parenth -= 1
                scan += 1
            the_sqrts.append( sexpr[pos:scan] )
    # 2/ Simplifications of the sqrt's
    new_expr = ""    # will contain the result
    pos0 = 0
    for i, pos in enumerate(pos_sqrts):
        # radcan is called on each sqrt:
        x = SR(the_sqrts[i])
        argum = x.operands()[0] # the argument of sqrt 
        den = argum.denominator()
        if den != 1:  # the argument of sqrt is a fraction
            num = argum.numerator()
            if num < 0 or den < 0:
                x = sqrt(-num) / sqrt(-den)  # new equivalent expression for x
        simpl = SR(x._maxima_().radcan())
        # the absolute value of radcan's output is taken, the call to simplify()
        # taking into account possible assumptions regarding the sign of simpl:
        ssimpl = str(abs(simpl).simplify())
        # search for abs(1/sqrt(...)) term to simplify it into 1/sqrt(...):
        pstart = ssimpl.find('abs(1/sqrt(')
        if pstart != -1:
            ssimpl = ssimpl[:pstart] + ssimpl[pstart+3:] # getting rid of 'abs'
        new_expr += sexpr[pos0:pos] + '(' + ssimpl + ')'
        pos0 = pos + len(the_sqrts[i])
    new_expr += sexpr[pos0:]
    return SR(new_expr)

def simplify_abs_trig(expr):
    r"""
    Simplify abs(sin(...)) in symbolic expressions.

    EXAMPLES::

        sage: forget()  # for doctests only
        sage: M = Manifold(3, 'M')
        sage: X.<x,y,z> = M.chart(r'x y:(0,pi) z:(-pi/3,0)')
        sage: X.coord_range()
        x: (-oo, +oo); y: (0, pi); z: (-1/3*pi, 0)

    Since ``x`` spans all `\RR`, no simplification of ``abs(sin(x))``
    occurs, while ``abs(sin(y))`` and ``abs(sin(3*z))`` are correctly
    simplified, given that `y \in (0,\pi)` and `z \in (-\pi/3,0)`::

        sage: from sage.geometry.manifolds.utilities import simplify_abs_trig
        sage: simplify_abs_trig( abs(sin(x)) + abs(sin(y)) + abs(sin(3*z)) )
        abs(sin(x)) + sin(y) - sin(3*z)

    Note that neither Sage's function
    :meth:`~sage.symbolic.expression.Expression.simplify_trig` nor
    :meth:`~sage.symbolic.expression.Expression.simplify_full`
    works in this case::

        sage: s = abs(sin(x)) + abs(sin(y)) + abs(sin(3*z))
        sage: s.simplify_trig()
        abs(4*cos(z)^2 - 1)*abs(sin(z)) + abs(sin(x)) + abs(sin(y))
        sage: s.simplify_full()
        abs(4*cos(z)^2 - 1)*abs(sin(z)) + abs(sin(x)) + abs(sin(y))

    despite the following assumptions hold::

        sage: assumptions()
        [x is real, y is real, y > 0, y < pi, z is real, z > -1/3*pi, z < 0]

    Additional checks are::

        sage: simplify_abs_trig( abs(sin(y/2)) )  # shall simplify
        sin(1/2*y)
        sage: simplify_abs_trig( abs(sin(2*y)) )  # must not simplify
        abs(sin(2*y))
        sage: simplify_abs_trig( abs(sin(z/2)) )  # shall simplify
        -sin(1/2*z)
        sage: simplify_abs_trig( abs(sin(4*z)) )  # must not simplify 
        abs(sin(4*z))

    """
    from sage.symbolic.ring import SR
    from sage.symbolic.constants import pi
    sexpr = str(expr)
    if 'abs(sin(' not in sexpr:  # nothing to simplify
        return expr
    tp = []
    val = []
    for pos in range(len(sexpr)):
        if sexpr[pos:pos+8] == 'abs(sin(':
            # finding the end of abs argument:
            scan = pos+4 # start of abs
            parenth = 1
            while parenth != 0:
                if sexpr[scan] == '(': parenth += 1
                if sexpr[scan] == ')': parenth -= 1
                scan += 1
            pos_abs_end = scan
            # finding the end of sin argument:
            scan = pos+8 # start of sin
            parenth = 1
            while parenth != 0:
                if sexpr[scan] == '(': parenth += 1
                if sexpr[scan] == ')': parenth -= 1
                scan += 1
            pos_sin_end = scan
            # if the abs contains only the sinus, the simplification can be tried:
            if pos_sin_end == pos_abs_end-1:
                tp.append(pos)
                val.append( sexpr[pos:pos_abs_end] )
    simp = []
    for v in val:
        # argument of the sinus:
        sx = v[8:-2]
        x = SR(sx)
        if x>=0 and x<=pi:
            simp.append('sin(' + sx + ')')
        elif x>=-pi and x<=0:
            simp.append('(-sin(' + sx + '))')
        else:
            simp.append(v)  # no simplification is applicable
    nexpr = ""
    pos0 = 0
    for i, pos in enumerate(tp):
        nexpr += sexpr[pos0:pos] + simp[i]
        pos0 = pos + len(val[i])
    nexpr += sexpr[pos0:]
    return SR(nexpr)



def simplify_chain(expr):
    r"""
    Apply a chain of simplications to a symbolic expression.

    This is the simplification chain used in calculus involving functions
    of coordinates in a given chart, as implemented in
    :class:`~sage.geometry.manifolds.chart.FunctionChart`.

    The chain is formed by the following functions, called
    successively:

    #. :meth:`~sage.symbolic.expression.Expression.simplify_factorial`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_trig`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_rational`
    #. :func:`simplify_sqrt_real`
    #. :func:`simplify_abs_trig`
    #. :meth:`~sage.symbolic.expression.Expression.canonicalize_radical`
       (for Sage >= 6.5) or
       :meth:`~sage.symbolic.expression.Expression.simplify_radical` (for
       Sage < 6.5)
    #. :meth:`~sage.symbolic.expression.Expression.simplify_log`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_rational`
    #. :meth:`~sage.symbolic.expression.Expression.simplify_trig`

    """
    expr = expr.simplify_factorial()
    expr = expr.simplify_trig()
    expr = expr.simplify_rational()
    expr = simplify_sqrt_real(expr)
    expr = simplify_abs_trig(expr)
    # In Sage 6.5, simplify_radical() has been renamed canonicalize_radical()
    #  (cf. http://trac.sagemath.org/11912): 
    if version == '6.5':
        expr = expr.canonicalize_radical()
    else:
        expr = expr.simplify_radical()
    expr = expr.simplify_log('one')
    expr = expr.simplify_rational()
    expr = expr.simplify_trig()
    return expr

def set_axes_labels(graph, xlabel, ylabel, zlabel, **kwds):
    r"""
    Set axes labels for a 3D graphics object.

    This is a workaround for the lack of axes labels in Sage 3D plots; it
    sets the labels as text3d objects at locations determined from the
    bounding box of the graphic object ``graph``.

    INPUT:

    - ``graph`` -- a 3D graphic object, as an instance of
      :class:`~sage.plot.plot3d.base.Graphics3d`
    - ``xlabel`` -- string for the x-axis label
    - ``ylabel`` -- string for the y-axis label
    - ``zlabel`` -- string for the z-axis label
    - ``**kwds`` -- options (e.g. color) for text3d

    OUTPUT:

    - the 3D graphic object with text3d labels added.

    """
    from sage.plot.plot3d.shapes2 import text3d
    xmin, ymin, zmin = graph.bounding_box()[0]
    xmax, ymax, zmax = graph.bounding_box()[1]
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin
    x1 = xmin + dx / 2
    y1 = ymin + dy / 2
    z1 = zmin + dz / 2
    xmin1 = xmin - dx / 20
    xmax1 = xmax + dx / 20
    ymin1 = ymin - dy / 20
    zmin1 = zmin - dz / 20
    graph += text3d('  ' + xlabel, (x1, ymin1, zmin1), **kwds)
    graph += text3d('  ' + ylabel, (xmax1, y1, zmin1), **kwds)
    graph += text3d('  ' + zlabel, (xmin1, ymin1, z1), **kwds)
    return graph
