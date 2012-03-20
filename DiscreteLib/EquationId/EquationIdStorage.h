
#pragma once

namespace DiscreteLib
{

/**
 * \brief Interface class for mapped address
 */
class IEquationIdStorage
{
public:
	virtual ~IEquationIdStorage() {};
    virtual bool hasKey(size_t pt_id) const = 0;
    virtual bool hasValue(size_t eqs_id) const = 0;
    virtual void key_range(size_t &i_start, size_t &i_end) const = 0;
    virtual void activate(size_t pt_id, bool b) = 0;
    virtual bool isActive(size_t pt_id) const = 0;
    virtual void set(size_t pt_id, long eqs_id) = 0;
    virtual size_t setAll(size_t address_start, size_t dn_pt=1) = 0;
    virtual size_t size() const = 0;
    virtual long address(size_t key_id) const = 0;
    virtual long key(long address_id) const = 0;
    virtual bool isSequential() const = 0;
};

} //end

