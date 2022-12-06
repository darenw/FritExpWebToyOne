// Functions to compute fractionally iterated exponentials, the Szekeres way
// 
// Daren Scot Wilson
// See PDF paper at https://github.com/darenw/FRITEXP/tree/master/doc
// or slides at
// https://www.slideshare.net/DarenWilson1/generalizing-addition-and-multiplication-to-an-operator-parametrized-by-a-real-number 



/*
   evalSzekeresPolys(u) returns the values of the Szekeres polynomials for a given 
   value of u.  These values are used as coefficients in a power series in x
   to compute g_u(x).
*/

const SzekPolys = [
   [ 0 ],
   [ 1 ],
   [ 0,   1/2  ],
   [ 0,  -1/12,   1/4  ],
   [ 0,   1/48,  -5/48,  1/8 ],
   [ 0,  -1/180,  1/24,  -13/144,  1/16  ],
   [ 0,  11/8640, -91/5760, 89/1728,  -77/1152,  1/32 ]
];


function poly(coefs, u)  {
    return coefs.reduceRight( (acc, coef) => coef+u*acc, 0);
}


function evalSzekeresPolys(u)  {
    return SzekPolys.map( sp =>  poly(sp, u) )
}


/* 
    Use this for -0.5 < u < 0.5,   -.1 < x < .1   (rough guidelines, cheating is allowed)
    Computes g_u(x) when both u and x are "small"
    g(x) = exp(x)-1  and g_u(x) is fractional iteration order u of g

    "Small" is undefined, but easy to get a feel for after performing some computations
    for any particular purpose.
    Recommendation: keep u in the
    range -0.5 to 0.5, and x to -0.1 to 0.1   

    Better (maybe): keep u in -0.25 to +0.25, repeat smallg() to effectively double u.

    Args:
        x (real scalar or np.array): main arg to feed the function
        u (real):   order of iteration of the function g(x) == exp(x)-1 
    Returns:
        real
    effectively double u.  Keep x in -0.5 to 0.5 or smaller.

*/

function gsmall(u, x)   {
    let xcoefs = evalSzekeresPolys(u);
    return poly(xcoefs, x);
}




const BIG_WONT_OVERFLOW = 709.7;  // for IEEE 64bit float, exp(this value) won't overflow
const SMALL_FOR_GUSMALL = 0.4;


function fritexp_core(xcoefs, u_int, x, small_limit=SMALL_FOR_GUSMALL) {
    
    let count =0;
    while (x < BIG_WONT_OVERFLOW)  {
        x = Math.exp(x);
        count +=1;
    }
    
    while (x> SMALL_FOR_GUSMALL)  {
        x = Math.log(1+x);
        count -= 1;
    }
    
    let y = poly(xcoefs, x)
    
    while (y < BIG_WONT_OVERFLOW)   {
        y = Math.exp(y) - 1;
        count += 1;
    }
    
    while (count - u_int > 0)  {
        if (y>0) {
            y = Math.log(y);
            count -= 1;
        } else {
            return NaN;
        }
    }
    
    return y;
}


/** 
 * Main function most users will want. 
 * Computes exp_u(x) for a given u and x.
 * Simple, trouble-free pure function.
 * May return useless values when it should return NAN, such as for log(-1)
 * @param  {Float} u     order of iteration
 * @param  {Float} x     value to give to the fractionally iterated exponential
 * @return  {Float}      value of exp_u(x)
 */
export function fritexp(u,x)    {
    return fritexp_function_factory(u)(x);
}



/**
 * This factory creates a function with some of the internal stuff precalculated.
 * Use that function to calculate faster many exp_u(x) for the same x 
 * @param  {Float}    The order of iteration, u 
 * @return {Function} Function that takes x, returns exp_u(x)
 */
export function fritexp_function_factory(u, small_limit=SMALL_FOR_GUSMALL)  {
    let u_int = Math.round(u); 
    let u_frac = u - u_int;
    let xcoefs = evalSzekeresPolys(u_frac);
    const f = function (x)  {
              return fritexp_core(xcoefs, u_int, x, small_limit);
    };
    return f;
}



//export  fritexp ;
//export { gsmall, fritexp, FritexpFunctionFactory }
