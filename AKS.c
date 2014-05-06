#include<flint/flint.h>
#include<flint/fmpz.h>
#include<flint/fmpz_mod_poly.h>
#include<flint/arith.h>
#include<math.h>
#include<stdio.h>
static __inline__ void MULTIPLICATIVE_ORDER(fmpz_t out,fmpz_t n,fmpz_t k){
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
    fmpz_t temp1;fmpz_init(temp1);
    /*Step1*/
    mpz_t k;mpz_init(k);fmpz_get_mpz(k,n);
    if(mpz_perfect_power_p(k)!=0){
            fmpz_clear(temp1);mpz_clear(k);
            return 0;
    }
    mpz_clear(k);
    /*Step2*/
    fmpz_t r;fmpz_init(r);
    ulong c;//rewrite c,change its type from fmpz to ulong
    c=powl(log2l(fmpz_get_d(n)),2);//c=[log(n)^2]
    for(fmpz_set_ui(r,2);;fmpz_add_ui(r,r,1)){
        MULTIPLICATIVE_ORDER(temp1,r,n);
        if(fmpz_cmp_ui(temp1,c)>0)break;
    }
    /*Step3*/
    fmpz_t a;fmpz_init(a);
    for(fmpz_one(a);fmpz_cmp(a,r)<=0;fmpz_add_ui(a,a,1)){
        fmpz_gcd(temp1,a,n);
        if(fmpz_cmp(temp1,n)<0&&fmpz_cmp_ui(temp1,1)>0){
            fmpz_clear(temp1);fmpz_clear(a);fmpz_clear(r);
            return 0;
        }
    }
    /*Step4*/
    if(fmpz_cmp(n,r)<=0){
            fmpz_clear(temp1);fmpz_clear(a);fmpz_clear(r);
            return 1;
    }
    /*Step5*/
    fmpz_clear(temp1);
    arith_euler_phi(temp1,r);
    c=ceil(sqrtl(fmpz_get_ui(temp1)))*log2l(fmpz_get_ui(n));//Recalculate c, thanks to Dana Jacobsen
    fmpz_mod_poly_t e1,e2;
    fmpz_mod_poly_init(e1,n);fmpz_mod_poly_init(e2,n);
    ulong n_ui=fmpz_get_ui(n);//Thanks to Denis Kryskov
    ulong r_ui=fmpz_get_ui(r);//Thanks to Denis Kryskov
    //move the definition of fmpz_mod_poly out of the for loop
    fmpz_mod_poly_t p     ;fmpz_mod_poly_init(p     ,n);
    fmpz_mod_poly_set_coeff_ui(p     ,1             , 1);
    fmpz_mod_poly_t q     ;fmpz_mod_poly_init(q     ,n);
    fmpz_mod_poly_set_coeff_ui(q     ,n_ui % r_ui, 1);//Thanks to Denis Kryskov,this modification really helps.
    fmpz_mod_poly_t modulo;fmpz_mod_poly_init(modulo,n);
    fmpz_mod_poly_set_coeff_ui(modulo,0             ,n_ui-1); //Thanks to Denis Kryskov
    printf("%llu\n",c);
    fmpz_mod_poly_set_coeff_ui(modulo,r_ui          , 1);//Thanks to Denis Kryskov
    for(fmpz_one(a);fmpz_cmp_ui(a,c)<=0;fmpz_add_ui(a,a,1)){
        fmpz_mod_poly_set_coeff_fmpz(p   ,0             , a);
        fmpz_mod_poly_set_coeff_fmpz(q   ,0             , a);//q=x^n+a or x^(n-r)+a
        fmpz_mod_poly_divrem(e1,e2,q,modulo);//e2=q mod modulo,e1 just takes the place here
        fmpz_mod_poly_powmod_fmpz_binexp(e1,p,n,modulo);//e1=p^n mod modulo
        if(fmpz_mod_poly_equal(e1,e2)==0){
                fmpz_clear(temp1);fmpz_clear(a);fmpz_clear(r);
                fmpz_mod_poly_clear(modulo);fmpz_mod_poly_clear(p);fmpz_mod_poly_clear(q);fmpz_mod_poly_clear(e1);fmpz_mod_poly_clear(e2);
                return 0;
        }
    }
    fmpz_mod_poly_clear(modulo);fmpz_mod_poly_clear(p);fmpz_mod_poly_clear(q);
    fmpz_clear(temp1);fmpz_clear(a);fmpz_clear(r);
    fmpz_mod_poly_clear(e1);fmpz_mod_poly_clear(e2);
    return 1;
}


int main(){
    fmpz_t test;fmpz_init(test);
    while(fmpz_read(test)){
        printf("%d\n",AKS(test));
    }
    fmpz_clear(test);
    system("pause");
    return 0;
}


