#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pymatgen as pmg
import numpy as np
import pycrystem as pc


# In[2]:


from pycrystem.utils.pyprismatic_io_utils import generate_pyprismatic_input as gpi


# In[3]:


Ga = pmg.Element("Ga")
As = pmg.Element("As")
lattice = pmg.Lattice.cubic(5.65)

structure = pmg.Structure.from_spacegroup("F-43m",lattice, [Ga,As], [[0, 0, 0],[0.5,0.5,0.5]]) #arg1


# In[4]:


gpi(structure,"test.xyz")


# In[4]:


structure

