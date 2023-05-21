#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "bn.h"
#include "bn_m.h"
#include "hcpy.h"

typedef long long Long;
#define BN_SIZE_T (BN_ARRAY_SIZE * sizeof(DTYPE))


UINT uint_pow(UINT x, UINT y){
	UINT i=0;
	UINT r=1;
	for(; i<y; i++) 
		r *= x;
	
	return r;
}

typedef struct k_v{
	// 引数 c 
	bn x_c, y_c, re_n, re2_n;
	
	// bn val 
	bn m, m2, a, b, c, d, ac, db, ad_bc, pow_t;	
	
	// uint val
	UINT m_n, m2_n, pow_m2;
	
	// bn k_mul val
	bn k_mul;

	// bn add val
	bn ab_add, cd_add, ac_db_add;

}k_v;

k_v	*K_VAL(void){
	k_v *v = (k_v*)malloc(sizeof(k_v));
	if(v == NULL){
		fprintf(stderr,"v malloc error");
		return NULL;
	}

	return v;
}

void karatuba_mul_b(bn *x, bn *y, bn *r){
	// 各種変数を生成
	k_v *k_val = K_VAL();
	memset(k_val, 0, sizeof(k_v));
	bn one, tow, third, teen; // int .. 0, 3, 10;
	
	// 変数を設定
	bn_from_int(&one,	1);
	bn_from_int(&tow,	2);
	bn_from_int(&third, 3);
	bn_from_int(&teen, 10);
	
	// 引数をコピー
	memcpy(&k_val->x_c, x, BN_SIZE_T);
	memcpy(&k_val->y_c, y, BN_SIZE_T);

	// 桁数を取得 
	UINT x_len = bn_dec_len(&k_val->x_c);
	UINT y_len = bn_dec_len(&k_val->y_c);
	
	printf("x_len or y_len : %d, %d\n", x_len, y_len);
	// 一桁である場合、引数をmul
	if(x_len < 1 || y_len < 1){
		bn_mul_s(x, y, r);
		free(k_val);
		return;
	}

	// 最大桁数を設定
	if(x_len < y_len) k_val->m_n = y_len;
	else		      k_val->m_n = x_len;

	k_val->m2_n = (k_val->m_n + 2) / 3 * 3 / 2;
	k_val->pow_m2 = uint_pow(10, k_val->m2_n);	
	
	bn_from_int(&k_val->m2, k_val->m2_n);
	bn_from_int(&k_val->pow_t, k_val->pow_m2);
	
	bn_div(&k_val->x_c, &k_val->pow_t, &k_val->a); // a = x / pow(10, m2_n);
	bn_mod(&k_val->x_c, &k_val->pow_t, &k_val->b); // b = x % pow(10, m2_n);
	bn_div(&k_val->y_c, &k_val->pow_t, &k_val->c); // c = y / pow(10, m2_n);	
	bn_mod(&k_val->y_c, &k_val->pow_t, &k_val->d); // d = y % pow(10, m2_n);
	
	#ifdef DEBUG
	int debug_a = bn_to_int(&k_val->a);
	int debug_b = bn_to_int(&k_val->b);
	int debug_c = bn_to_int(&k_val->c);
	int debug_d = bn_to_int(&k_val->d);
	
	printf("\n----\na : %d\nb : %d\nc : %d\nd : %d\n----\n",
		debug_a,debug_b,debug_c,debug_d);
	#endif 
	
	karatuba_mul_b(&k_val->a, &k_val->c, &k_val->ac); // ac = a * c
	karatuba_mul_b(&k_val->d, &k_val->b, &k_val->db); // db = d * b;

	bn_add(&k_val->a, &k_val->b, &k_val->ab_add); // ab_add = a + c;
	bn_add(&k_val->c, &k_val->d, &k_val->cd_add); // cd_add = c + d;
	
	#ifdef DEBUG
	puts("ab_add or cd_add");
	bn_from_dec_print(&k_val->ab_add);
	bn_from_dec_print(&k_val->cd_add);
	#endif 

	karatuba_mul_b(&k_val->ab_add, &k_val->cd_add, &k_val->k_mul); // k_mul = ab_add * cd_add;
	
	bn_sub(&k_val->k_mul, &k_val->ac, &k_val->k_mul); // k_mul = k_mul - ac;
	bn_sub(&k_val->k_mul, &k_val->db, &k_val->k_mul); // k_mul = k_mul - db;
	
	memcpy(&k_val->ad_bc, &k_val->k_mul, BN_SIZE_T); // ad_bc = k_mul;

	#ifdef DEBUG 
	int debug_ac = bn_to_int(&k_val->ac);
	int debug_db = bn_to_int(&k_val->db);
	int debug_ad_bc = bn_to_int(&k_val->ad_bc);

	printf("\n----\nac : %d\ndb : %d\nad_bc : %d\n----\n",
			debug_ac,debug_db,debug_ad_bc);
	#endif	
	
	karatuba_mul_b(&k_val->m2, &tow, &k_val->re_n);			// re_n = m2 * 2;
	bn_pow_s(&teen, &k_val->re_n, &k_val->re_n);			// re_n = pow(10, re_n);
	bn_mul_s(&k_val->ac, &k_val->re_n, &k_val->re_n);		// re_n = ac * re_n;
	
	bn_pow_s(&teen, &k_val->m2, &k_val->re2_n);				// re2_n = pow(10, m2);
	bn_mul_s(&k_val->ad_bc, &k_val->re2_n, &k_val->re2_n);	// re2_n = ad_bc * re2_n;
	bn_add(&k_val->db, &k_val->re2_n, &k_val->re2_n);		// re2_n = db + re2_n
	
	bn_add(&k_val->re_n, &k_val->re2_n, r);	// r = re_n + re2_n;
	free(k_val);
	return;
}


void karatuba_mul(bn *x, bn *y, bn *r){
	// 各種変数を生成
	k_v *k_val = K_VAL();
	memset(k_val, 0, sizeof(k_v));
	bn one, tow, third, teen; // int .. 0, 3, 10;
	
	// 変数を設定
	bn_from_int(&one,	1);
	bn_from_int(&tow,	2);
	bn_from_int(&third, 3);
	bn_from_int(&teen, 10);
	
	// 引数をコピー
	memcpy(&k_val->x_c, x, BN_SIZE_T);
	memcpy(&k_val->y_c, y, BN_SIZE_T);
	
	
	// 桁数を取得 
	UINT x_len = bn_dec_len(&k_val->x_c);
	UINT y_len = bn_dec_len(&k_val->y_c);
	
	//printf("x_len or y_len : %d, %d\n", x_len, y_len);
	// 一桁である場合、引数をmul
	if(x_len <= 1 || y_len <= 1){
		bn_mul_s(x, y, r);
		free(k_val);
		return;
	}

	// 最大桁数を設定
	if(x_len < y_len) k_val->m_n = y_len;
	else		      k_val->m_n = x_len;

	k_val->m2_n = (k_val->m_n + 2) / 3 * 3 / 2;
	k_val->pow_m2 = uint_pow(10, k_val->m2_n);	
	
	bn_from_int(&k_val->m2, k_val->m2_n);
	bn_from_int(&k_val->pow_t, k_val->pow_m2);
	bn_div(&k_val->x_c, &k_val->pow_t, &k_val->a); // a = x / pow(10, m2_n);
	bn_mod(&k_val->x_c, &k_val->pow_t, &k_val->b); // b = x % pow(10, m2_n);
	bn_div(&k_val->y_c, &k_val->pow_t, &k_val->c); // c = y / pow(10, m2_n);	
	bn_mod(&k_val->y_c, &k_val->pow_t, &k_val->d); // d = y % pow(10, m2_n);
	
	#ifdef DEBUG
	int debug_a = bn_to_int(&k_val->a);
	int debug_b = bn_to_int(&k_val->b);
	int debug_c = bn_to_int(&k_val->c);
	int debug_d = bn_to_int(&k_val->d);
	
	printf("\n----\na : %d\nb : %d\nc : %d\nd : %d\n----\n",
		debug_a,debug_b,debug_c,debug_d);
	#endif 
	
	karatuba_mul(&k_val->a, &k_val->c, &k_val->ac); // ac = a * c
	karatuba_mul(&k_val->d, &k_val->b, &k_val->db); // db = d * b;

	bn_add(&k_val->a, &k_val->b, &k_val->ab_add); // ab_add = a + c;
	bn_add(&k_val->c, &k_val->d, &k_val->cd_add); // cd_add = c + d;
	
	#ifdef DEBUG
	puts("ab_add or cd_add");
	bn_from_dec_print(&k_val->ab_add);
	bn_from_dec_print(&k_val->cd_add);
	#endif 

	karatuba_mul(&k_val->ab_add, &k_val->cd_add, &k_val->k_mul); // k_mul = ab_add * cd_add;
	
	bn_sub(&k_val->k_mul, &k_val->ac, &k_val->k_mul); // k_mul = k_mul - ac;
	bn_sub(&k_val->k_mul, &k_val->db, &k_val->k_mul); // k_mul = k_mul - db;
	
	memcpy(&k_val->ad_bc, &k_val->k_mul, BN_SIZE_T); // ad_bc = k_mul;

	#ifdef DEBUG 
	int debug_ac = bn_to_int(&k_val->ac);
	int debug_db = bn_to_int(&k_val->db);
	int debug_ad_bc = bn_to_int(&k_val->ad_bc);

	printf("\n----\nac : %d\ndb : %d\nad_bc : %d\n----\n",
			debug_ac,debug_db,debug_ad_bc);
	#endif	
	
	bn_mul_s(&k_val->m2, &tow, &k_val->re_n);		// re_n = m2 * 2;
	bn_pow_s(&teen, &k_val->re_n, &k_val->re_n);		// re_n = pow(10, re_n);
	bn_mul_s(&k_val->ac, &k_val->re_n, &k_val->re_n);	// re_n = ac * re_n;
	
	bn_pow_s(&teen, &k_val->m2, &k_val->re2_n);				// re2_n = pow(10, m2);
	bn_mul_s(&k_val->ad_bc, &k_val->re2_n, &k_val->re2_n);	// re2_n = ad_bc * re2_n;
	bn_add(&k_val->db, &k_val->re2_n, &k_val->re2_n);		// re2_n = db + re2_n
	
	bn_add(&k_val->re_n, &k_val->re2_n, r);	// r = re_n + re2_n;
	free(k_val);
	return;
}
