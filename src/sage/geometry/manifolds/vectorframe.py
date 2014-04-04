r"""
Vector Frames

The class :class:`VectorFrame` implements vector frames on differentiable 
manifolds over `\RR`. By 'vector frame' it is meant a field on a manifold M that
provides, at each point p in M, a vector basis of the tangent space at p. 

A derived class of :class:`VectorFrame` is :class:`CoordFrame`; it regards the 
vector frames associated with a chart, i.e. the so-called coordinate bases. 

The vector frame duals, i.e. the coframes, are implemented via the class
:class:`CoFrame`. The derived class :class:`CoordCoFrame` is devoted to 
coframes deriving from a chart. 

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2013, 2014): initial version

EXAMPLES:
    
    Setting a vector frame on a 3-dimensional manifold::
    
        sage: m = Manifold(3, 'M')
        sage: c_xyz = m.chart('x y z')
        sage: e = VectorFrame(m, 'e') ; e
        vector frame (M, (e_0,e_1,e_2))
        sage: latex(e)
        \left(M ,\left(e_0,e_1,e_2\right)\right)

    The first frame defined on a manifold is its default frame; in the present
    case it is the coordinate frame defined when introducing the chart c_xyz::
    
        sage: m.default_frame()
        coordinate frame (M, (d/dx,d/dy,d/dz))
        
    The default frame can be changed via the method
    :meth:`Manifold.set_default_frame`::
    
        sage: m.set_default_frame(e)
        sage: m.default_frame()
        vector frame (M, (e_0,e_1,e_2))

    The elements of a vector frame are vector fields on the manifold::
    
        sage: e.vec
        (vector field 'e_0' on the 3-dimensional manifold 'M', vector field 'e_1' on the 3-dimensional manifold 'M', vector field 'e_2' on the 3-dimensional manifold 'M')   
        
    Each element can be accessed by its index::
    
        sage: e[0]
        vector field 'e_0' on the 3-dimensional manifold 'M'
            
    The index range depends on the starting index defined on the manifold::

        sage: m = Manifold(3, 'M', start_index=1)
        sage: c_xyz = m.chart('x y z')
        sage: e = VectorFrame(m, 'e')              
        sage: e.vec
        (vector field 'e_1' on the 3-dimensional manifold 'M', vector field 'e_2' on the 3-dimensional manifold 'M', vector field 'e_3' on the 3-dimensional manifold 'M')
        sage: e[1], e[2], e[3]
        (vector field 'e_1' on the 3-dimensional manifold 'M', vector field 'e_2' on the 3-dimensional manifold 'M', vector field 'e_3' on the 3-dimensional manifold 'M')
    
    Let us check that the vector fields e(i) are the frame vectors from
    their components w.r.t. to the frame e::
    
        sage: e[1].comp(e)[:]
        [1, 0, 0]
        sage: e[2].comp(e)[:]
        [0, 1, 0]
        sage: e[3].comp(e)[:]
        [0, 0, 1]
    
    Defining a vector frame on a manifold automatically creates the dual 
    coframe, which bares the same name (here e)::
    
        sage: m.coframes
        [coordinate coframe (M, (dx,dy,dz)), coframe (M, (e^1,e^2,e^3))]
        sage: f = m.coframes[1] ; f
        coframe (M, (e^1,e^2,e^3))
   
    Each element of the coframe is a 1-form::
   
        sage: f[1], f[2], f[3]
        (1-form 'e^1' on the 3-dimensional manifold 'M',
        1-form 'e^2' on the 3-dimensional manifold 'M',
        1-form 'e^3' on the 3-dimensional manifold 'M')
        sage: latex(f[1]), latex(f[2]), latex(f[3])
        (e^1, e^2, e^3)

    Let us check that the coframe (e^i) is indeed the dual of the vector 
    frame (e_i)::
    
        sage: f[1](e[1]) # the 1-form e^1 applied to the vector field e_1
        scalar field 'e^1(e_1)' on the 3-dimensional manifold 'M'
        sage: f[1](e[1]).expr() # the explicit expression of e^1(e_1)
        1
        sage: f[1](e[1]).expr(), f[1](e[2]).expr(), f[1](e[3]).expr()
        (1, 0, 0)
        sage: f[2](e[1]).expr(), f[2](e[2]).expr(), f[2](e[3]).expr()
        (0, 1, 0)
        sage: f[3](e[1]).expr(), f[3](e[2]).expr(), f[3](e[3]).expr()
        (0, 0, 1)
    
    The coordinate frame associated to spherical coordinates of the 
    sphere `S^2`::
    
        sage: m = Manifold(2, 'S^2', start_index=1)
        sage: c_spher.<th,ph> = m.chart(r'th:[0,pi]:\theta ph:[0,2*pi):\phi')
        sage: b = m.default_frame() ; b
        coordinate frame (S^2, (d/dth,d/dph))
        sage: b[1]
        vector field 'd/dth' on the 2-dimensional manifold 'S^2'
        sage: b[2]
        vector field 'd/dph' on the 2-dimensional manifold 'S^2'

    The orthonormal frame constructed from the coordinate frame::
    
        sage: change_frame = AutomorphismField(m)
        sage: change_frame[:] = [[1,0], [0, 1/sin(th)]]
        sage: e = b.new_frame(change_frame, 'e') ; e
        vector frame (S^2, (e_1,e_2))
        sage: e[1][:]
        [1, 0]
        sage: e[2][:]
        [0, 1/sin(th)]
        
    The change-of-frame matrices::
    
        sage: m.frame_change(c_spher.frame, e)        
        field of tangent-space automorphisms on the 2-dimensional manifold 'S^2'
        sage: m.frame_change(c_spher.frame, e)[:]
        [        1         0]
        [        0 1/sin(th)]
        sage: m.frame_change(e, c_spher.frame)
        field of tangent-space automorphisms on the 2-dimensional manifold 'S^2'
        sage: m.frame_change(e, c_spher.frame)[:]
        [      1       0]
        [      0 sin(th)]

"""

#******************************************************************************
#       Copyright (C) 2013, 2014 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2013, 2014 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.tensor.modules.free_module_basis import FreeModuleBasis, FreeModuleCoBasis

class VectorFrame(FreeModuleBasis):
    r"""
    Class for vector frames on a differentiable manifold over `\RR`. 
    
    By *vector frame*, it is meant a field `e` on some open domain `U` of a 
    manifold `S` endowed with a mapping `\Phi: U\rightarrow V` to a 
    parallelizable domain `V` of a manifold `M` such that for each `p\in U`, 
    `e(p)` is a vector basis of the tangent space `T_{\Phi(p)}M`
    
    The standard case of a vector frame *on* `U` corresponds to `S=M`, `U=V`
    and `\Phi = \mathrm{Id}`. 
    
    For each instanciation of a vector frame, a coframe is automatically 
    created, as an instance of the class :class:`CoFrame`. 
    
    INPUT:
    
    - ``domain`` -- manifold domain on which the vector frame is defined
    - ``symbol`` -- (default: None) a letter (of a few letters) to denote a 
      generic vector of the frame; can be set to None if the parameter
      ``from_frame`` is filled.
    - ``latex_symbol`` -- (default: None) symbol to denote a generic vector of
      the frame; if None, the value of ``symbol`` is used. 
    - ``dest_map`` -- (default: None) destination map `\Phi:\ U \rightarrow V` 
      (type: :class:`~sage.geometry.manifolds.diffmapping.DiffMapping`); 
      if none is provided, the identity is assumed (case of a vector frame *on* 
      `U`)
    - ``from_frame`` -- (default: None) vector frame `\tilde e` on the codomain 
      `V` of the destination map `\Phi`; the frame `e` = ``self`` is then 
      constructed so that `\forall p \in U, e(p) = \tilde e(\Phi(p))`


    EXAMPLES:

    Setting a vector frame on a 3-dimensional manifold::
    
        sage: m = Manifold(3, 'M')
        sage: c_xyz = m.chart('x y z')
        sage: e = VectorFrame(m, 'e')
        sage: e
        vector frame (M, (e_0,e_1,e_2))
        sage: latex(e)
        \left(M ,\left(e_0,e_1,e_2\right)\right)

    The LaTeX symbol can be specified::
    
        sage: e = VectorFrame(m, 'E', r"\epsilon")
        sage: latex(e)
        \left(M ,\left(\epsilon_0,\epsilon_1,\epsilon_2\right)\right)

    
    """
    def __init__(self, domain, symbol=None, latex_symbol=None, dest_map=None,
                 from_frame=None):
        self.domain = domain
        self.dest_map = dest_map
        self.from_frame = from_frame
        self.manifold = domain.manifold
        if symbol is None:
            if from_frame is None:
                raise TypeError("Some frame symbol must be provided.") 
            symbol = 'X'  # provisory symbol
        FreeModuleBasis.__init__(self, domain.vector_field_module(dest_map), 
                                 symbol, latex_symbol=latex_symbol)
        # Redefinition of the name and the LaTeX name:
        if from_frame is None:
            self.name = "(" + self.domain.name + ", " + self.name + ")"
            self.latex_name = r"\left(" + self.domain.latex_name + ", " + \
                          self.latex_name + r"\right)"
        else:
            if not from_frame.domain.is_subdomain(dest_map.codomain):
                raise ValueError("The domain of the frame 'from_frame' is " + 
                                 "not included in the codomain of the " + 
                                 "destination map.")
            n = self.fmodule.rank()
            for i in range(n):
                self.vec[i].name = from_frame.vec[i].name
                self.vec[i].latex_name = from_frame.vec[i].latex_name
            self.name = "(" + self.domain.name + ", (" + \
                        ",".join([self.vec[i].name for i in range(n)]) + "))"
            self.latex_name = r"\left(" + self.domain.latex_name + \
                        r" ,\left(" + \
                        ",".join([self.vec[i].latex_name for i in range(n)])+ \
                        r"\right)\right)"
            self.symbol = from_frame.symbol
            self.latex_symbol = from_frame.latex_symbol
        # The frame is added to the domain's set of frames, as well as to all 
        # the superdomains' sets of frames; moreover the first defined frame 
        # is considered as the default one
        for sd in self.domain.superdomains:
            for other in sd.frames:
                if repr(self) == repr(other):
                    raise ValueError("The " + str(self) + " already exist on" +
                                     " the " + str(sd))
            sd.frames.append(self)
            if sd.def_frame is None: 
                sd.def_frame = self
        if dest_map is None:
            # The frame is added to the list of the domain's covering frames:
            self.domain.covering_frames.append(self)
        #
        # Dual coframe 
        self.coframe = self.dual_basis()  # self.coframe = a shortcut for self._dual_basis
        #
        # Derived quantities:
        self._structure_coef = None
        # Initialization of the set of frames that are restrictions of the
        # current frame to subdomains of the frame domain:
        self.subframes = set([self]) 
        # Initialization of the set of frames which the current frame is a 
        # restriction of:
        self.superframes = set([self]) 


    ###### Methods that must be redefined by derived classes of FreeModuleBasis ######

    def _repr_(self):
        r"""
        String representation of the object.
        """
        description = "vector frame " + self.name
        if self.dest_map is not None:
            description += " with values on the " + str(self.dest_map.codomain)
        return description
        

    def _init_dual_basis(self):
        r""" 
        Construct the basis dual to ``self``.
        
        OUTPUT:
        
        - instance of :class:`CoFrame` representing the dual of
          ``self``
        
        """
        return CoFrame(self, self.symbol, latex_symbol=self.latex_symbol)

    ###### End of methods redefined by derived classes ######

        
    def new_frame(self, change_of_frame, symbol, latex_symbol=None):
        r"""
        Define a new vector frame from the current one. 
        
        The new vector frame is defined on the same domain as ``self`` from
        a field of automorphisms. 
        
        INPUT:
        
        - ``change_of_frame`` -- instance of 
          :class:`~sage.geometry.rank2field.AutomorphismField`
          describing the automorphism `P` that relates the current frame 
          `(e_i)` (described by ``self``) to the new frame `(n_i)` according
          to `n_i = P(e_i)`
        - ``symbol`` -- a letter (of a few letters) to denote a generic vector of
          the frame
        - ``latex_symbol`` -- (default: None) symbol to denote a generic vector of
          the frame; if None, the value of ``symbol`` is used. 
          
        OUTPUT:
        
        - the new frame `(n_i)`, as an instance of :class:`VectorFrame`
        
        EXAMPLES:
        
        Frame resulting from a pi/3-rotation in the Euclidean plane::
        
            sage: m = Manifold(2,'R^2')
            sage: c_xy = m.chart('x y')
            sage: e = VectorFrame(m, 'e') ; m.set_default_frame(e)
            sage: m.frame_changes
            {}
            sage: rot = AutomorphismField(m)
            sage: rot[:] = [[sqrt(3)/2, -1/2], [1/2, sqrt(3)/2]]
            sage: n = e.new_frame(rot, 'n')
            sage: n[0][:]
            [1/2*sqrt(3), 1/2]
            sage: n[1][:]
            [-1/2, 1/2*sqrt(3)]
            sage: a =  m.frame_change(e,n)
            sage: a[:]
            [1/2*sqrt(3)        -1/2]
            [        1/2 1/2*sqrt(3)]
            sage: a == rot
            True
            sage: a is rot
            False
            sage: a.components
            {vector frame (R^2, (e_0,e_1)): 2-indices components w.r.t. the vector frame (R^2, (e_0,e_1)),
            vector frame (R^2, (n_0,n_1)): 2-indices components w.r.t. the vector frame (R^2, (n_0,n_1))}
            sage: a.comp(n)[:]
            [1/2*sqrt(3)        -1/2]
            [        1/2 1/2*sqrt(3)]
            sage: a1 = m.frame_change(n,e)
            sage: a1[:]
            [1/2*sqrt(3)         1/2]
            [       -1/2 1/2*sqrt(3)]
            sage: a1 == rot.inverse()
            True
            sage: a1 is rot.inverse()
            False
            sage: e[0].comp(n)[:]
            [1/2*sqrt(3), -1/2]
            sage: e[1].comp(n)[:]
            [1/2, 1/2*sqrt(3)]
  
        """
        return self.new_basis(change_of_frame, symbol, latex_symbol=latex_symbol)
        
    def new_subframe(self, domain, symbol, latex_symbol=None):
        r"""
        Construct a subframe.
        
        If ``self`` is a vector frame defined on the domain U, a subframe
        is the restriction of ``self`` to a subdomain V of U.
        
        INPUT:
        
        - ``domain`` -- subdomain `V` of the current frame domain `U` 
        - ``symbol`` -- a letter (of a few letters) to denote a generic vector of
          the frame
        - ``latex_symbol`` -- (default: None) symbol to denote a generic vector of
          the frame; if None, the value of ``symbol`` is used. 
        
        OUTPUT:
        
        - the subframe, as an instance of :class:`VectorFrame`. 

        """
        if not domain.is_subdomain(self.domain):
            raise TypeError("The argument 'domain' must be a subdomain of " + 
                            " the frame domain.")
        #!# to be changed: des_map should be the restriction of self.dest_map 
        # to V:
        res = VectorFrame(domain, symbol, latex_symbol=latex_symbol, 
                          dest_map=self.dest_map)
        # Update of superframes and subframes:
        res.superframes.update(self.superframes)
        for sframe in self.superframes:
            sframe.subframes.add(res)
        return res
    
    def structure_coef(self):
        r"""
        Evaluate the structure coefficients associated to the vector frame. 
        
        `n` being the manifold's dimension, the structure coefficients of the
        vector frame `(e_i)` are the `n^3` scalar fields `C^k_{\ \, ij}` 
        defined by 
        
        .. MATH::
            
            [e_i, e_j] = C^k_{\ \, ij} e_k
            
        OUPUT:
        
        - the structure coefficients `C^k_{\ \, ij}`, as an instance of 
          :class:`CompWithSym` with 3 indices ordered as `(k,i,j)`. 
          
        EXAMPLE:
        
        Structure coefficients of the orthonormal frame associated to
        spherical coordinates in the Euclidean space `R^3`::
        
            sage: m = Manifold(3, 'R^3', '\RR^3', start_index=1)
            sage: c_spher.<r,th,ph> = m.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi')
            sage: ch_frame = AutomorphismField(m) 
            sage: ch_frame[1,1], ch_frame[2,2], ch_frame[3,3] = 1, 1/r, 1/(r*sin(th))
            sage: m.frames
            [coordinate frame (R^3, (d/dr,d/dth,d/dph))]
            sage: e = c_spher.frame.new_frame(ch_frame, 'e')
            sage: e[1][:]  # components of e_1 in the manifold's default frame (d/dr, d/dth, d/dth)
            [1, 0, 0]
            sage: e[2][:]
            [0, 1/r, 0]
            sage: e[3][:]
            [0, 0, 1/(r*sin(th))]
            sage: c = e.structure_coef() ; c
            3-indices components w.r.t. the vector frame (R^3, (e_1,e_2,e_3)), with antisymmetry on the index positions (1, 2)
            sage: c[:]
            [[[0, 0, 0], [0, 0, 0], [0, 0, 0]],
             [[0, -1/r, 0], [1/r, 0, 0], [0, 0, 0]],
             [[0, 0, -1/r], [0, 0, -cos(th)/(r*sin(th))], [1/r, cos(th)/(r*sin(th)), 0]]]
            sage: c[2,1,2]  # C^2_{12}
            -1/r
            sage: c[3,1,3]  # C^3_{13}
            -1/r
            sage: c[3,2,3]  # C^3_{23}
            -cos(th)/(r*sin(th))


        """
        from component import CompWithSym
        if self._structure_coef is None:
            self._structure_coef = CompWithSym(self, 3, antisym=(1,2))
            si = self.manifold.sindex
            nsi = si + self.manifold.dim
            for k in range(si,nsi):
                ce_k = self.coframe.form[k-si]
                for i in range(si, nsi):
                    e_i = self.vec[i-si]
                    for j in range(i+1, nsi):
                        e_j = self.vec[j-si]
                        self._structure_coef[[k,i,j]] = ce_k(e_j.lie_der(e_i))
        return self._structure_coef
            
#******************************************************************************

class CoordFrame(VectorFrame):
    r"""
    Class for coordinate frames on a differentiable manifold over `\RR`. 
    
    By 'coordinate frame', it is meant a vector frame on a manifold M that 
    is associated to a coordinate system (chart) on M. 
    
    INPUT:
    
    - ``chart`` -- the chart defining the coordinates

    EXAMPLES:

    The coordinate frame associated to spherical coordinates of the 
    sphere `S^2`::
    
        sage: m = Manifold(2, 'S^2', start_index=1)
        sage: m.chart(r'th:[0,pi]:\theta ph:[0,2*pi):\phi')
        chart (S^2, (th, ph))
        sage: b = m.default_frame()
        sage: b
        coordinate frame (S^2, (d/dth,d/dph))
        sage: b[1]
        vector field 'd/dth' on the 2-dimensional manifold 'S^2'
        sage: b[2]
        vector field 'd/dph' on the 2-dimensional manifold 'S^2'
        sage: latex(b)
        \left(S^2 ,\left(\frac{\partial}{\partial \theta },\frac{\partial}{\partial \phi }\right)\right)
 
    """
    def __init__(self, chart):
        from sage.misc.latex import latex
        from chart import Chart
        if not isinstance(chart, Chart):
            raise TypeError("The first argument must be a chart.")
        self.chart = chart
        VectorFrame.__init__(self, chart.domain, symbol='X') # 'X' = provisory symbol
        n = self.manifold.dim
        for i in range(n):
            self.vec[i].name = "d/d" + str(self.chart.xx[i])
            self.vec[i].latex_name = r"\frac{\partial}{\partial" + \
                                     latex(self.chart.xx[i]) + r"}"
        self.name = "(" + self.domain.name + ", (" + \
                    ",".join([self.vec[i].name for i in range(n)]) + "))"
        self.latex_name = r"\left(" + self.domain.latex_name + r" ,\left(" + \
                       ",".join([self.vec[i].latex_name for i in range(n)])+ \
                       r"\right)\right)"
        self.symbol = self.name
        self.latex_symbol = self.latex_name


    ###### Methods that must be redefined by derived classes of FreeModuleBasis ######

    def _repr_(self):
        r"""
        String representation of the object.
        """
        return "coordinate frame " + self.name

    def _init_dual_basis(self):
        r""" 
        Construct the basis dual to ``self``.
        
        OUTPUT:
        
        - instance of :class:`FreeModuleCoBasis` representing the dual of
          ``self``
        
        """
        return CoordCoFrame(self)

    ###### End of methods redefined by derived classes ######

    def structure_coef(self):
        r"""
        Returns the structure coefficients associated to the vector frame. 
        
        `n` being the manifold's dimension, the structure coefficients of the
        vector frame `(e_i)` are the `n^3` scalar fields `C^k_{\ \, ij}` 
        defined by 
        
        .. MATH::
            
            [e_i, e_j] = C^k_{\ \, ij} e_k
        
        In the present case, where `(e_i)` is a coordinate frame, 
        `C^k_{\ \, ij}=0`. 
        
        OUPUT:
        
        - the structure coefficients `C^k_{\ \, ij}`, as a vanishing instance 
          of :class:`CompWithSym` with 3 indices ordered as `(k,i,j)`
          
        EXAMPLE:
        
        Structure coefficients of the coordinate frame associated to
        spherical coordinates in the Euclidean space `R^3`::
        
            sage: m = Manifold(3, 'R^3', r'\RR^3', start_index=1)
            sage: c_spher = m.chart(r'r:[0,+oo) th:[0,pi]:\theta ph:[0,2*pi):\phi')
            sage: b = m.default_frame() ; b
            coordinate frame (R^3, (d/dr,d/dth,d/dph))
            sage: c = b.structure_coef() ; c
            3-indices components w.r.t. the coordinate frame (R^3, (d/dr,d/dth,d/dph)), with antisymmetry on the index positions (1, 2)
            sage: c == 0
            True
        
        """
        from component import CompWithSym
        if self._structure_coef is None:
            self._structure_coef = CompWithSym(self, 3, antisym=(1,2))
            # A just created CompWithSym is zero
        return self._structure_coef
        

#******************************************************************************

class CoFrame(FreeModuleCoBasis):
    r"""
    Class for coframes on a differentiable manifold over `\RR`. 
    
    By 'coframe', it is meant a n-tuple of 1-forms on a manifold M that 
    provides, at each point p in M, a basis of the space dual to the tangent 
    space at p. 
    
    INPUT:
    
    - ``frame`` -- the vector frame dual to the coframe
    - ``symbol`` -- a letter (of a few letters) to denote a generic 1-form in
      the coframe
    - ``latex_symbol`` -- (default: None) symbol to denote a generic 1-form in
      the coframe; if None, the value of ``symbol`` is used. 

    EXAMPLES:

    Coframe on a 3-dimensional manifold::
    
        sage: m = Manifold(3, 'M', start_index=1)
        sage: c_xyz = m.chart('x y z')
        sage: v = VectorFrame(m, 'v')
        sage: e = CoFrame(v, 'e') ; e
        coframe (M, (e^1,e^2,e^3))

    The 1-forms composing the coframe are obtained via the () operator::
        
        sage: e[1], e[2], e[3]
        (1-form 'e^1' on the 3-dimensional manifold 'M',
         1-form 'e^2' on the 3-dimensional manifold 'M',
         1-form 'e^3' on the 3-dimensional manifold 'M')

    Checking that e is the dual of v::
    
        sage: e[1](v[1]).expr(), e[1](v[2]).expr(), e[1](v[3]).expr()
        (1, 0, 0)
        sage: e[2](v[1]).expr(), e[2](v[2]).expr(), e[2](v[3]).expr()
        (0, 1, 0)
        sage: e[3](v[1]).expr(), e[3](v[2]).expr(), e[3](v[3]).expr()
        (0, 0, 1)

    """
    def __init__(self, frame, symbol, latex_symbol=None):
        self.domain = frame.domain
        self.manifold = self.domain.manifold
        FreeModuleCoBasis.__init__(self, frame, symbol, 
                                   latex_symbol=latex_symbol)
        # Redefinition of the name and the LaTeX name:
        self.name = "(" + self.domain.name + ", " + self.name + ")"
        self.latex_name = r"\left(" + self.domain.latex_name + ", " + \
                          self.latex_name + r"\right)"
        # The coframe is added to the domain's set of coframes, as well as to 
        # all the superdomains' sets of coframes
        for sd in self.domain.superdomains:
            for other in sd.coframes:
                if repr(self) == repr(other):
                    raise ValueError("The " + str(self) + " already exist on" +
                                     " the " + str(sd))
            sd.coframes.append(self)
        
    
    def _repr_(self):
        r"""
        String representation of the object.
        """
        return "coframe " + self.name


#******************************************************************************

class CoordCoFrame(CoFrame):
    r"""
    Class for coordinate coframes on a differentiable manifold over `\RR`. 
    
    By 'coordinate coframe', it is meant the n-tuple of the differentials of 
    the coordinates belonging to a chart on a manifold.
    
    INPUT:
    
    - ``coord_frame`` -- coordinate frame dual to the coordinate coframe

    EXAMPLES:
    
    Coordinate coframe on a 3-dimensional manifold::
    
        sage: m = Manifold(3, 'M', start_index=1)
        sage: c_xyz = m.chart('x y z')
        sage: m.frames
        [coordinate frame (M, (d/dx,d/dy,d/dz))] 
        sage: m.coframes
        [coordinate coframe (M, (dx,dy,dz))]
        sage: dX = m.coframes[0] ; dX
        coordinate coframe (M, (dx,dy,dz))

    The 1-forms composing the coframe are obtained via the () operator::
    
        sage: dX[1]
        1-form 'dx' on the 3-dimensional manifold 'M'
        sage: dX[2]
        1-form 'dy' on the 3-dimensional manifold 'M'
        sage: dX[3]
        1-form 'dz' on the 3-dimensional manifold 'M'
        sage: dX[1][:]
        [1, 0, 0]
        sage: dX[2][:]
        [0, 1, 0]
        sage: dX[3][:]
        [0, 0, 1]

    The coframe is the dual of the coordinate frame::
    
        sage: e = c_xyz.frame ; e
        coordinate frame (M, (d/dx,d/dy,d/dz))
        sage: dX[1](e[1]).expr(), dX[1](e[2]).expr(), dX[1](e[3]).expr()
        (1, 0, 0)
        sage: dX[2](e[1]).expr(), dX[2](e[2]).expr(), dX[2](e[3]).expr()
        (0, 1, 0)
        sage: dX[3](e[1]).expr(), dX[3](e[2]).expr(), dX[3](e[3]).expr()
        (0, 0, 1)
    
    Each 1-form of a coordinate coframe is closed::
    
        sage: dX[1].exterior_der()
        2-form 'ddx' on the 3-dimensional manifold 'M'
        sage: dX[1].exterior_der() == 0
        True

    """
    def __init__(self, coord_frame):
        from sage.misc.latex import latex
        if not isinstance(coord_frame, CoordFrame):
            raise TypeError("The first argument must be a coordinate frame.")
        CoFrame.__init__(self, coord_frame, 'X') # 'X' = provisory symbol
        self.chart = coord_frame.chart
        n = self.manifold.dim
        for i in range(n):
            self.form[i].name = "d" + str(self.chart.xx[i])
            self.form[i].latex_name = r"\mathrm{d}" + latex(self.chart.xx[i])
        self.name = "(" + self.domain.name + ", (" + \
                    ",".join([self.form[i].name for i in range(n)]) + "))"
        self.latex_name = r"\left(" + self.domain.latex_name + \
                          r" ,\left(" + \
                ",".join([self.form[i].latex_name for i in range(n)])+ \
                          r"\right)\right)"

    def _repr_(self):
        r"""
        String representation of the object.
        """
        return "coordinate coframe " + self.name 
