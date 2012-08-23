#ifndef CHEMCOMP_H
#define CHEMCOMP_H

#include "chemconst.h"

namespace ogsChem
{

// decalre ChemPhase class here to avoid cross inclusion.
class ChemPhase;
class ChemReaction;

class ChemComp
{
public:
    // constructor and destructor
    /**
      * constructor of the class
      */
    ChemComp(void);
    /**
      * distructor of the class
      */
    ~ChemComp(void);

    /**
      * getter for mole mass in gram/mol.
      */
    double get_mole_mass(void) const;

    /**
      * getter for mole volume in m3.
      */
    double get_mole_volume(void) const;

    /**
      * getter for charge in the unit of equivalent.
      */
    double get_charge(void) const;

    double get_neg_gamma_1(void)
    {return neg_gamma_1;}

    double get_neg_gamma_2(void)
    {return neg_gamma_2;}

    /**
      * getter for chemical component name.
      */
    std::string& get_name(void)
    {return _name;}

    void set_name(std::string new_name)
    {_name = new_name;}

    void set_composition(std::string new_composition)
    {_composition = new_composition;}

    std::string get_composition(void)
    {return _composition;}

    void set_compTyps(Comp_Type new_type)
    {_mCompType = new_type;}

    void set_neg_gamma(double n_gamma_1, double n_gamma_2)
    {neg_gamma_1 = n_gamma_1; neg_gamma_2 = n_gamma_2;}

    void set_pReaction(ChemReaction* new_Reaction)
    {_mReaction = new_Reaction;}

    ChemReaction* get_pReaction(void)
    {return _mReaction;}

private:
    /**
      * The name of the chemical component.
      */
    std::string _name;

    /**
      * The molecular composition of the component.
      */
    std::string _composition;

    /**
      * Molar mass of the component, unit in g/mol.
      */
    double _mole_mass;

    /**
      * Molar volume of the component, unit in m3/mol.
      */
    double _mole_volume;

    /**
      * charge of the component, unit in eq/mol,
      * e.g. calsium cation Ca2+ has a charge of +2.0 eq/mol
      */
    double _charge;

    /**
      * Gibbs formation energy of the component.
      */
    double _G0;

    /**
      * pointer of phase this component belongs to.
      */
    ChemPhase* _mPhase;

    /**
      * pointer of phase this component belongs to.
      */
    ChemReaction* _mReaction;

    /**
      * type of phase this component belongs to.
      */
    Comp_Type _mCompType;

    /**
      * negative gamma values
      */
    double neg_gamma_1, neg_gamma_2;

};

} // end of namespace

#endif // CHEMCOMP_H
