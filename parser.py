#!/usr/bin/env python3

"""
A parser of SMILES chemical notation using pyparsing module the EBNF of SMILES
is taken from:
https://metamolecular.com/cheminformatics/smiles/railroad-diagram/
Code Refined by NKUCodingCat@Nov8,2019
"""

import pyparsing as pp

# Grammar definition
isotope = pp.Word(pp.nums) #
atomclass = pp.Regex(':[0-9]+') #
bond = pp.oneOf(['-','=','#','$',':','/','\\','.'])
organicsymbol = pp.oneOf(['B','Br','C','Cl','N','O','P','S','F','I'])
aromaticsymbol = pp.oneOf(['b','c','n','o','p','s'])
elementsymbol = pp.oneOf(['H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne', 'Na', 'Mg', 
                          'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca', 'Sc', 'Ti', 'V',  'Cr',
                          'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
                          'Rb', 'Sr', 'Y',  'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                          'In', 'Sn', 'Sb', 'Te',  'I', 'Xe', 'Cs', 'Ba', 'Hf', 'Ta', 'W',  'Re',
                          'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr',
                          'Ra', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Fl', 'Lv',
                          'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er',
                          'Tm', 'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'U',  'Np', 'Pu', 'Am', 'Cm', 'Bk',
                          'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'])
hcount = pp.Literal('H') + pp.Optional(pp.Regex('[0-9]')) # 
ringclosure = pp.Optional( pp.Literal('%') + pp.Regex('[1-9]')) + pp.Regex('[1-9]')
charge = (pp.Literal('-') + pp.Optional( pp.oneOf( ['-', ] + list(map(str, range(0, 17))) ) ) ) ^\
         (pp.Literal('+') + pp.Optional( pp.oneOf( ['+', ] + list(map(str, range(0, 17))) ) ) ) #
chiralclass = pp.Optional(
                 pp.Literal('@') ^ pp.Literal('@@') ^ \
                (pp.Literal('@') + ( \
                  ( pp.oneOf(['TH', 'AL'])  + pp.Regex('[1-2]') ) ^\
                  ( pp.Literal('SP') + pp.Regex('[1-3]') ) ^\
                  ( pp.Literal('TB') + pp.oneOf(list(map(str, range(1, 21)))) ) ^\
                  ( pp.Literal('OH') + pp.oneOf(list(map(str, range(1, 31)))) ) 
                )) 
            ) #
atomspec = pp.Literal('[') +\
                pp.Optional(isotope) + \
                ( pp.Literal('se') ^ pp.Literal('as') ^ aromaticsymbol ^ elementsymbol ^ pp.Literal('*') ) + \
                pp.Optional(chiralclass)+ pp.Optional(hcount) + pp.Optional(charge) + pp.Optional(atomclass) + \
           pp.Literal(']') #
atom = organicsymbol ^ aromaticsymbol ^ pp.Literal('*') ^ atomspec #
chain = pp.OneOrMore(pp.Optional(bond) + ( atom ^ ringclosure )) #
## This looks fucked up
smiles = pp.Forward()
branch = pp.Forward()
smiles << atom + pp.ZeroOrMore(chain ^ branch)
branch << pp.Literal('(') + pp.Optional(bond) +  pp.OneOrMore(smiles) + pp.Literal(')')



def IsValidSMILES(text):
    """
    A simple SMILES validator
    """
    is_valid = False
    results = smiles.parseString(text, parseAll=True)
    if results:
        is_valid = True
        return(is_valid)
    return is_valid

if __name__ == '__main__':
    astr = 'CC(CC)CO[Na+]'
    print(IsValidSMILES(astr)) 
    astr = 'F[C@H](C1=C([NH3+])C=C(F)C=C1)C2=CC([O-])=CC=C2'
    print(IsValidSMILES(astr)) 
    astr = 'F[C@H](C1=CC=C(F)C=C1)C2=CC=CC=C2'
    print(IsValidSMILES(astr)) 
    astr = 'F[14C@@H](C1=CC=C(F)C=C1)C2=CC=CC=C2'
    print(IsValidSMILES(astr)) 
    
