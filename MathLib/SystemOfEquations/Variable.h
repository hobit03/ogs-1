
#pragma once

#include <string>

namespace MathLib
{

struct Variable
{
	size_t id;
	std::string name;
	size_t n_dof;

	Variable(size_t i, size_t n) : id(i), n_dof(n) {};
};


} //end
