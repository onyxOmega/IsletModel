#include "BCell.h"

#include <iostream>

void BCell::initialize()
{	
	testDouble = 1000.0;
}

double BCell::getTestDouble() const
{
	return testDouble;
}
