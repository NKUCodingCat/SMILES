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
elementsymbol = pp.oneOf(['Al','Am','Sb','Ar','33','At','Ba','Bk','Be','Bi',
                           'Bh','B','Br','Cd','Ca','Cf','C','Ce','Cs','Cl','Cr',
                           'Co','Cu','Cm','Ds','Db','Dy','Es','Er','Eu','Fm',
                           'F','Fr','Gd','Ga','Ge','Au','Hf','Hs','He','Ho','H',
                           'In','I','Ir','Fe','Kr','La','Lr','Pb','Li','Lu','Mg',
                           'Mn','Mt','Md','Hg','Mo','Nd','Ne','Np','Ni','Nb','N',
                           'No','Os','O','Pd','P','Pt','Pu','Po','K','Pr','Pm',
                           'Pa','Ra','Rn','Re','Rh','Rg','Rb','Ru','Rf','Sm',
                           'Sc','Sg','Se','Si','Ag','Na','Sr','S','Ta','Tc',
                           'Te','Tb','Tl','Th','Tm','Sn','Ti','W','Uub','Uuh',
                           'Uuo','Uup','Uuq','Uus','Uut','Uuu','U','V','Xe','Yb',
                           'Y','Zn','Zr',])
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
    
