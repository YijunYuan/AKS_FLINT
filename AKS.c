#include<flint/flint.h>
#include<flint/fmpz.h>
#include<flint/fmpz_mod_poly.h>
#include<flint/arith.h>
#include<stdio.h>
#include<math.h>
void MULTIPLICATIVE_ORDER(fmpz_t out,fmpz_t n,fmpz_t k){
    fmpz_gcd(out,n,k);
    if(!fmpz_equal_ui(out,1)){
        fmpz_set_si(out,-1);return;
    }
    fmpz_t i;fmpz_init(i);
    for(fmpz_one(i);;fmpz_add_ui(i,i,1)){
        fmpz_powm(out,k,i,n);
        if(fmpz_equal_ui(out,1)){
                fmpz_set(out,i);
                fmpz_clear(i);
                return;
        }
    }
}
int AKS(fmpz_t n){
    if(fmpz_equal_ui(n,2)==1)return 1;
    if(fmpz_divisible_si(n,2)==1)return 0;
    fmpz_t temp1;fmpz_init(temp1);
	if (fmpz_equal_ui(n, 1) == 1){
            fmpz_clear(temp1);
            return 0;
    }
	/*Step1*/
    mpz_t k;mpz_init(k);fmpz_get_mpz(k,n);
	if(mpz_perfect_power_p(k)!=0){
            fmpz_clear(temp1);mpz_clear(k);
            return 0;
    }
	mpz_clear(k);
	/*Step2*/
	fmpz_t c,r;fmpz_init(c);fmpz_init(r);
	fmpz_set_d(c,sqrt(fmpz_dlog(n)));//c=[sqrt(log(n))]
    for(fmpz_set_ui(r,2);;fmpz_add_ui(r,r,1)){
        MULTIPLICATIVE_ORDER(temp1,r,n);
        if(fmpz_cmp(temp1,c)>0)break;
    }
    /*Step3*/
    fmpz_t a;fmpz_init(a);
    for(fmpz_one(a);fmpz_cmp(a,r)<=0;fmpz_add_ui(a,a,1)){
        fmpz_gcd(temp1,a,n);
        if(fmpz_cmp(temp1,n)<0&&fmpz_cmp_ui(temp1,1)>0){
            fmpz_clear(temp1);fmpz_clear(a);fmpz_clear(c);fmpz_clear(r);
            return 0;
        }
    }
    /*Step4*/
    if(fmpz_cmp(n,r)<=0){
            fmpz_clear(temp1);fmpz_clear(a);fmpz_clear(c);fmpz_clear(r);
            return 1;
    }
    /*Step5*/
    fmpz_set_d(c,fmpz_dlog(n));//c=log(n)
    arith_euler_phi(temp1,r);//temp1=euler(r)
    fmpz_mul(c,c,temp1);//c=log(n)*euler(r)
    fmpz_clear(temp1);
    fmpz_mod_poly_t e;fmpz_mod_poly_init(e,n);
    for(fmpz_one(a);fmpz_cmp(a,c)<=0;fmpz_add_ui(a,a,1)){
        fmpz_mod_poly_t modulo;fmpz_mod_poly_init(modulo,n);
        fmpz_mod_poly_t p     ;fmpz_mod_poly_init(p     ,n);
        fmpz_mod_poly_t q     ;fmpz_mod_poly_init(q     ,n);
        fmpz_mod_poly_set_coeff_ui(modulo,0             ,-1);
        fmpz_mod_poly_set_coeff_ui(modulo,fmpz_get_ui(r), 1);//modulo=x^r-1
        fmpz_mod_poly_set_coeff_ui(p     ,1             , 1);
        fmpz_mod_poly_set_coeff_fmpz(p   ,0             , a);//p=x+a
        fmpz_mod_poly_set_coeff_ui(q     ,fmpz_get_ui(n), 1);
        fmpz_mod_poly_set_coeff_fmpz(q   ,0             , a);//q=x^n+a
        fmpz_mod_poly_powmod_fmpz_binexp(p,p,n,modulo);//p=p^n mod modulo
        fmpz_mod_poly_divrem(e,q,q,modulo);//q=q mod modulo
        if(fmpz_mod_poly_equal(p,q)==0){
                fmpz_clear(temp1);fmpz_clear(a);fmpz_clear(c);fmpz_clear(r);
                fmpz_mod_poly_clear(modulo);fmpz_mod_poly_clear(p);fmpz_mod_poly_clear(q);fmpz_mod_poly_clear(e);
                return 0;
        }
        fmpz_mod_poly_clear(modulo);fmpz_mod_poly_clear(p);fmpz_mod_poly_clear(q);
    }
    fmpz_clear(temp1);fmpz_clear(a);fmpz_clear(c);fmpz_clear(r);
    fmpz_mod_poly_clear(e);
	return 1;
}
int main(){
    fmpz_t test;fmpz_init(test);
    while(1){
        fmpz_read(test);
        printf("%d\n",AKS(test));
    }
    fmpz_clear(test);
    return 0;
}
