r"""
Tensors in index notation

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version

"""
#******************************************************************************
#       Copyright (C) 2014 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2014 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.sage_object import SageObject

class TensorWithIndices(SageObject):
    r"""
    Tensors in index notation
    """
    def __init__(self, tensor, indices):
        self._tensor = tensor
        self._changed = False # indicates whether self contains an altered 
                              # version of the original tensor (True if 
                              # symmetries are indicated in the indices)
        # Suppress all '{' and '}' comming from LaTeX notations:
        indices = indices.replace('{','').replace('}','')
        # Suppress the first '^':
        if indices[0] == '^':
            indices = indices[1:]
        if '^' in indices:
            raise TypeError("The contravariant indices must be placed first.")
        con_cov = indices.split('_')
        con = con_cov[0]
        # Contravariant indices
        # ---------------------
        #  Search for (anti)symmetries:
        if '(' in con:
            sym1 = con.index('(')
            sym2 = con.index(')')-2
            if con.find('(', sym1+1) != -1 or '[' in con:
                raise NotImplementedError("Multiple symmetries are not " + 
                                          "treated yet.")
            self._tensor = self._tensor.symmetrize(range(sym1, sym2+1))
            self._changed = True # self does no longer contain the original tensor
            con = con.replace('(','').replace(')','')
        if '[' in con:
            sym1 = con.index('[')
            sym2 = con.index(']')-2
            if con.find('[', sym1+1) != -1 or '(' in con:
                raise NotImplementedError("Multiple symmetries are not " + 
                                          "treated yet.")
            self._tensor = self._tensor.antisymmetrize(range(sym1, sym2+1))
            self._changed = True # self does no longer contain the original tensor
            con = con.replace('[','').replace(']','')
        if len(con) != self._tensor._tensor_type[0]:
            raise TypeError("Number of contravariant indices not compatible " + 
                            "with the tensor type.")
        self._con = con
        # Covariant indices
        # -----------------
        if len(con_cov) == 1:
            if tensor._tensor_type[1] != 0:
                raise TypeError("Number of covariant indices not compatible " + 
                                "with the tensor type.")
            self._cov = ''
        elif len(con_cov) == 2:
            cov = con_cov[1]
            #  Search for (anti)symmetries:
            if '(' in cov:
                sym1 = cov.index('(')
                sym2 = cov.index(')')-2
                if cov.find('(', sym1+1) != -1 or '[' in cov:
                    raise NotImplementedError("Multiple symmetries are not " + 
                                              "treated yet.")
                csym1 = sym1 + self._tensor._tensor_type[0]
                csym2 = sym2 + self._tensor._tensor_type[0]
                self._tensor = self._tensor.symmetrize(range(csym1, csym2+1))
                self._changed = True # self does no longer contain the original 
                                     # tensor
                cov = cov.replace('(','').replace(')','')
            if '[' in cov:
                sym1 = cov.index('[')
                sym2 = cov.index(']')-2
                if cov.find('[', sym1+1) != -1 or '(' in cov:
                    raise NotImplementedError("Multiple symmetries are not " + 
                                              "treated yet.")
                csym1 = sym1 + self._tensor._tensor_type[0]
                csym2 = sym2 + self._tensor._tensor_type[0]
                self._tensor = self._tensor.antisymmetrize(range(csym1, csym2+1))
                self._changed = True # self does no longer contain the original 
                                     # tensor
                cov = cov.replace('[','').replace(']','')
            if len(cov) != tensor._tensor_type[1]:
                raise TypeError("Number of covariant indices not compatible " + 
                                "with the tensor type.")
            self._cov = cov
        else:
            raise TypeError("Two many '_' in the list of indices.")

    def _repr_(self):
        r"""
        String representation of the object.
        """
        if self._tensor._name is not None:
            name = self._tensor._name
        else:
            name = 'X'
        if self._con == '':
            return name + '_' + self._cov
        elif self._cov == '':
            return name + '^' + self._con
        else:
            return name + '^' + self._con + '_' + self._cov

    def update(self):
        r"""
        Return the tensor contains in ``self`` if it differs from that used
        for creating ``self``, otherwise return ``self``. 
        """
        if self._changed:
            return self._tensor
        else:
            return self

    def __mul__(self, other):
        r"""
        Tensor contraction on specified indices
        """
        if not isinstance(other, TensorWithIndices):
            raise TypeError("The second item of * must be a tensor with " + 
                            "specified indices.")
        contraction_pairs = []
        for ind in self._con:
            if ind != '.':
                if  ind in other._cov:
                    pos1 = self._con.index(ind)
                    pos2 = other._tensor._tensor_type[0] + other._cov.index(ind)
                    contraction_pairs.append((pos1, pos2))
                if  ind in other._con:
                    raise TypeError("The index " + str(ind) + " appears twice "
                                    + "in a contravariant position.")
        for ind in self._cov:
            if ind != '.':
                if ind in other._con:
                    pos1 = self._tensor._tensor_type[0] + self._cov.index(ind)
                    pos2 = other._con.index(ind)
                    contraction_pairs.append((pos1, pos2))
                if ind in other._cov:
                    raise TypeError("The index " + str(ind) + " appears twice "
                                    + "in a covariant position.")
        if contraction_pairs == []:
            # No contraction is performed: the tensor product is returned
            return self._tensor * other._tensor
        if len(contraction_pairs) > 1:
            raise NotImplementedError("Multiple contractions are not " + 
                                      "implemented yet.")
        pos1 = contraction_pairs[0][0]
        pos2 = contraction_pairs[0][1]
        return self._tensor.contract(pos1, other._tensor, pos2)

        
