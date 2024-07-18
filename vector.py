from numpy import array, zeros, dot

class Vector:
    def __init__(self, n):
        self.vec = zeros((n,), dtype=int)
        self.sq_length = 0
    
    def __len__(self):
        return self.sq_length**0.5
    
    def __getitem__(self, index):
        return self.vec[index]
    
    def __setitem__(self, index, value):
        temp = self.vec[index]
        self.vec[index] = value
        # squared norm
        self.sq_length += value*value - temp*temp
    
    def setarray(self, tab):
        self.vec = tab
        self.sq_length = dot(tab, tab)
        
    def setlist(self, lst):
        assert len(self.vec) == len(lst)
        norm = 0
        for i in range(len(self.vec)):
            self.vec[i] = lst[i]
            norm += lst[i]*lst[i]
        self.sq_length = norm
            
        
    def __eq__(self, other):
        assert len(self.vec) == len(other.vec)
        n = len(self.vec)
        i = 0
        j = 0
        while i < n and self.vec[i] == other.vec[i]:
            i += 1
        while j < n and self.vec[j] == -other.vec[j]:
            j += 1
        return (i == n) or (j == n)
    
    def __lt__(self, other):
        assert len(self.vec) == len(other.vec)
        return self.sq_length < other.sq_length
    
    def __gt__(self, other):
        assert len(self.vec) == len(other.vec)
        return self.sq_length > other.sq_length
    
    def __mul__(self, other):
        assert len(self.vec) == len(other.vec)
        return dot(self.vec, other.vec)

    def __sub__(self, other):
        """
        compute the shorter vector between u+v and u-v
        """
        assert len(self.vec) == len(other.vec)
        n = len(self.vec)
        scalar = self.__mul__(other)
        norm = 0
        v_new = Vector(n)
        if scalar < 0:
            for i in range(n):
                v_new.vec[i] = self.vec[i] + other.vec[i]
                norm += v_new.vec[i]*v_new.vec[i]
            v_new.sq_length = norm
        else:
            for i in range(n):
                v_new.vec[i] = self.vec[i] - other.vec[i]
                norm += v_new.vec[i]*v_new.vec[i]
            v_new.sq_length = norm
        return v_new
    
    def __repr__(self):
        return f"vector : {self.vec}\n norm² : {self.sq_length}"