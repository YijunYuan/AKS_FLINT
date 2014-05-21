AKS_FLINT
=========

Implementation of AKS Algorithm with FLINT

William Hart's FLINT is necessary. It's available at http://www.flintlib.org/

In theory, the max effective value of this implementation is the max value of long long unsigned int, but the real case still needs testing.

I'm glad to see that Denis Kryskov modified this code heavily and applied it into his Python library Razin. It's here https://github.com/krisk0/razin/blob/master/flint.binding/C/fmpz/AKS_trunc.c May be his version is better and more reliable.

Dana Jacobsen helped me to fix several serious BUGs, and FLINT-devel members also help a lot.
