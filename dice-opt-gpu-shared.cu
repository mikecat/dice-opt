#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <ctype.h>
#include <errno.h>

#define THREADS_ALIGN 32

/* strが有効なuint32_t型の範囲の正の整数かを返す */
int is_valid_pnum(const char* str) {
	char* end;
	long res;
	if (str == NULL) return 0;
	errno = 0;
	res = strtol(str, &end, 10);
	return *end == '\0' && errno == 0 && 1 <= res && (unsigned long)res <= UINT32_MAX;
}

/* オーバーフローに気をつけながら、size_tにuint32_tを掛ける */
size_t multiply_size(size_t size, uint32_t m) {
	if (size > SIZE_MAX / m) {
		fprintf(stderr, "size to allocate too big!\n");
		exit(1);
	}
	return size * m;
}

/* dst += src; 最後の繰り上がりの数を返す */
uint32_t add_nums(uint32_t* dst, const uint32_t* src, uint32_t length) {
	uint32_t carry = 0;
	uint32_t i;
	for (i = 0; i < length; i++) {
		uint32_t next_carry = (dst[i] > UINT32_MAX - src[i]);
		dst[i] += src[i];
		next_carry += (dst[i] > UINT32_MAX - carry);
		dst[i] += carry;
		carry = next_carry;
	}
	return carry;
}
__inline__ __device__ uint32_t add_nums_gpu(uint32_t* dst, const uint32_t* src, uint32_t length) {
	uint32_t carry = 0;
	uint32_t i;
	for (i = 0; i < length; i++) {
		uint32_t next_carry = (dst[i] > UINT32_MAX - src[i]);
		dst[i] += src[i];
		next_carry += (dst[i] > UINT32_MAX - carry);
		dst[i] += carry;
		carry = next_carry;
	}
	return carry;
}

/* out = in * mul; 最後の繰り上がりの数を返す */
uint32_t multiply(uint32_t* out, const uint32_t* in, uint32_t mul, uint32_t length) {
	uint32_t i;
	uint32_t carry = 0;
	for (i = 0; i < length; i++) {
		uint64_t next = (uint64_t)in[i] * mul + carry;
		out[i] = (uint32_t)(next & UINT32_C(0xffffffff));
		carry = (uint32_t)(next >> 32);
	}
	return carry;
}

/* 指定されたサイコロのパターン数を表すのに必要十分なDWORD数を返す */
/* 返す数はalignmentの倍数とする */
uint32_t get_required_dwords(
uint32_t dice_num, uint32_t dice_max_output, uint32_t alignment) {
	uint32_t* calc_buffer = (uint32_t*)malloc(sizeof(uint32_t));
	uint32_t calc_buffer_len = 1;
	uint32_t i;
	uint32_t carry;
	uint64_t ret;
	if (calc_buffer == NULL) {
		perror("get_required_bytes : malloc");
		exit(1);
	}
	calc_buffer[0] = 1;
	for (i = 0; i < dice_num; i++) {
		carry = multiply(calc_buffer, calc_buffer, dice_max_output, calc_buffer_len);
		if (carry > 0) {
			uint32_t* next_buffer = (uint32_t*)realloc(calc_buffer,
				sizeof(uint32_t) * ++calc_buffer_len);
			if (next_buffer == NULL) {
				perror("get_required_bytes : realloc");
				free(calc_buffer);
				exit(1);
			}
			calc_buffer = next_buffer;
			calc_buffer[calc_buffer_len - 1] = carry;
		}
	}
	free(calc_buffer);

	ret = (((uint64_t)calc_buffer_len + alignment - 1) / alignment)
		* alignment;
	if (ret > UINT32_MAX) {
		fprintf(stderr, "get_required_bytes : size too big\n");
		exit(1);
	}
	return (uint32_t)ret;
}

struct result_data {
	uint32_t dice_num, dice_max_output;
	uint32_t *puttern_count, *all_puttern_count;
};

uint32_t cmp_size, cmp_size2;
uint32_t *cmp_buffer, *cmp_buffer_a, *cmp_buffer_b, *cmp_buffer_temp;
int result_data_cmp(const void* x, const void* y) {
	const struct result_data *a = (const struct result_data*)x;
	const struct result_data *b = (const struct result_data*)y;
	uint32_t i;
	/* cmp_buffer_a = a->puttern_count * b->all_puttern_count */
	cmp_buffer_a[cmp_size] = multiply(cmp_buffer_a, a->puttern_count,
		b->all_puttern_count[0], cmp_size);
	for (i = 1; i < cmp_size; i++) {
		cmp_buffer_a[(size_t)cmp_size + i] = multiply(cmp_buffer_temp,
			a->puttern_count, b->all_puttern_count[i], cmp_size);
		cmp_buffer_a[(size_t)cmp_size + i] += add_nums(cmp_buffer_a + i,
			cmp_buffer_temp, cmp_size);
	}
	/* cmp_buffer_b = b->puttern_count * a->all_puttern_count */
	cmp_buffer_b[cmp_size] = multiply(cmp_buffer_b, b->puttern_count,
		a->all_puttern_count[0], cmp_size);
	for (i = 1; i < cmp_size; i++) {
		cmp_buffer_b[(size_t)cmp_size + i] = multiply(cmp_buffer_temp,
			b->puttern_count, a->all_puttern_count[i], cmp_size);
		cmp_buffer_b[(size_t)cmp_size + i] += add_nums(cmp_buffer_b + i,
			cmp_buffer_temp, cmp_size);
	}\
	/* 掛け算の結果の比較 (降順) */
	for (i = 0; i < cmp_size2; i++) {
		if (cmp_buffer_a[cmp_size2 - i - 1] > cmp_buffer_b[cmp_size2 - i - 1])
			return -1;
		if (cmp_buffer_a[cmp_size2 - i - 1] < cmp_buffer_b[cmp_size2 - i - 1])
			return 1;
	}
	/* サイコロを振る数の比較 (昇順) */
	if (a->dice_num > b->dice_num) return 1;
	if (a->dice_num < b->dice_num) return -1;
	/* サイコロで出る目の最大値の比較 (昇順) */
	if (a->dice_max_output > b->dice_max_output) return 1;
	if (a->dice_max_output< b->dice_max_output) return -1;
	return 0;
}

/* 多倍長整数をdouble型に変換する */
double num_to_double(const uint32_t* num, uint32_t length) {
	double res = 0;
	uint32_t i;
	for (i = 0; i < length; i++) {
		res = res * (UINT32_MAX + 1.0) + num[length - i - 1];
	}
	return res;
}

/* 多倍長整数を10進数で出力する。入力の多倍長整数は破壊される。 */
void print_num_and_destroy(uint32_t* num, uint32_t length) {
	static const uint32_t divisor = UINT32_C(1000000000);
	uint32_t remainder = 0;
	uint32_t i;
	int nonzero_exists = 0;
	for (i = 0; i < length; i++) {
		uint64_t current = ((uint64_t)remainder << 32) + num[length - i - 1];
		if ((num[length - i - 1] = (uint32_t)(current /divisor)) != 0) nonzero_exists = 1;
		remainder = (uint32_t)(current % divisor);
	}
	if (nonzero_exists) {
		print_num_and_destroy(num, length);
		printf("%09"PRIu32, remainder);
	} else {
		printf("%"PRIu32, remainder);
	}
}

/* 初期化 */
__global__ void init_dp(uint32_t* calculate_src, uint32_t required_dwords) {
	uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
	for (uint32_t j = 0; j < required_dwords; j++) {
		calculate_src[i * (size_t)required_dwords + j] = 0;
	}
	/* サイコロを0個振った時、合計は0の1通りのみ */
	if (i == 0) {
		calculate_src[0] = 1;
	}
}

__global__ void run_dp(uint32_t* calculate_src, uint32_t* calculate_dst,
uint32_t required_dwords, uint32_t dice_max) {
	uint32_t i = blockIdx.x * blockDim.x + threadIdx.x;
	/* その数について、パターン数を計算する */
	extern __shared__ uint32_t dst_buffer[];
	/* 足し算のために初期化する */
	for (uint32_t j = 0; j < required_dwords; j++) {
		dst_buffer[threadIdx.x * (size_t)required_dwords + j] = 0;
	}
	/* パターン数の計算のため、前のパターン数を足す */
	for (uint32_t j = 0; j <= dice_max; j++) {
		if (i > j) {
			add_nums_gpu(
				&dst_buffer[threadIdx.x * (size_t)required_dwords],
				&calculate_src[(i - j - 1) * (size_t)required_dwords],
				required_dwords);
		}
	}
	/* 結果を書き込む */
	for (uint32_t j = 0; j < required_dwords; j++) {
		calculate_dst[i * (size_t)required_dwords + j] =
			dst_buffer[threadIdx.x * (size_t)required_dwords + j];
	}
}

int main(int argc, char* argv[]) {
	/* 入力パラメータ */
	uint32_t dice_max_num, dice_max_output, target_value, output_num;
	/* パターン数を表すのに使うDWORD数 */
	uint32_t required_dwords;
	/* 計算バッファ */
	uint32_t *calculate_buffer, *calculate_src, *calculate_dst, *tmp;
	size_t calculate_buffer_size;
	size_t calculate_buffer_offset;
	/* 計算に用いる変数 */
	uint32_t dice_max_sum;
	uint32_t i, j, k;
	/* 結果の保存に用いる変数 */
	struct result_data *results, *current_result;
	uint32_t *result_putterns, *result_all_putterns;
	uint32_t result_count = 0;
	uint32_t *previous_all_putterns = NULL;

	/* 入力を読み取る */
	if ((argc != 4 && argc != 5) || !is_valid_pnum(argv[1]) ||
	!is_valid_pnum(argv[2]) || !is_valid_pnum(argv[3]) ||
	(argc >= 5 && !is_valid_pnum(argv[4]))) {
		fprintf(stderr,
			"Usage: %s max-num-of-dice max-output-of-a-die target-value [output-num]\n",
			argc >= 1 ? argv[0] : "dice-opt");
		return 1;
	}
	dice_max_num = (uint32_t)atol(argv[1]);
	dice_max_output = (uint32_t)atol(argv[2]);
	target_value = (uint32_t)atol(argv[3]);
	output_num = argc >= 5 ? (uint32_t)atol(argv[4]) : UINT32_C(10);
	printf("calculating putterns that becomes %"PRIu32"\n"
		"with max %"PRIu32" dice whose output are upto max %"PRIu32"\n",
		target_value, dice_max_num, dice_max_output);

	/* 計算バッファを確保する */
	if (dice_max_num > UINT32_MAX / dice_max_output) {
		fputs("parameter too big!\n", stderr);
		return 1;
	}
	dice_max_sum = dice_max_output * dice_max_num;
	if (dice_max_sum >= UINT32_MAX - ((THREADS_ALIGN - dice_max_sum % THREADS_ALIGN) % THREADS_ALIGN)) {
		fputs("parameter too big!\n", stderr);
		return 1;
	}
	if (target_value > dice_max_sum) {
		puts("there is no chance because target is too big");
		return 0;
	}
	dice_max_sum++;
	if (dice_max_sum % THREADS_ALIGN != 0) {
		uint32_t delta = THREADS_ALIGN - dice_max_sum % THREADS_ALIGN;
		if (dice_max_sum > UINT32_MAX - delta) {
			fputs("parameter too big!\n", stderr);
			return 1;
		}
		dice_max_sum += delta;
	}
	required_dwords = get_required_dwords(dice_max_num, dice_max_output, 1);
	calculate_buffer_offset = multiply_size(
			multiply_size(sizeof(uint32_t), required_dwords), dice_max_sum);
	calculate_buffer_size = multiply_size(calculate_buffer_offset, 2);
	cudaMalloc(&calculate_buffer, calculate_buffer_size);
	calculate_src = calculate_buffer;
	calculate_dst = calculate_buffer +
		(calculate_buffer_offset / sizeof(uint32_t));

	/* 結果の保存に用いるバッファの確保を行う */
	results = (struct result_data*)malloc(
		multiply_size(sizeof(struct result_data), dice_max_sum));
	if (results == NULL) {
		perror("failed to allocate result buffer");
		free(calculate_buffer);
		return 1;
	}
	result_putterns = (uint32_t*)malloc(
		multiply_size(
			multiply_size(sizeof(uint32_t), dice_max_sum), required_dwords));
	if (result_putterns == NULL) {
		perror("failed to allocate result puttern buffer");
		free(calculate_buffer);
		free(results);
		return 1;
	}
	result_all_putterns = (uint32_t*)malloc(
		multiply_size(
			multiply_size(sizeof(uint32_t), dice_max_sum), required_dwords));
	if (result_putterns == NULL) {
		perror("failed to allocate result all puttern buffer");
		free(calculate_buffer);
		free(results);
		free(result_putterns);
		return 1;
	}

	/* 比較用のメモリを確保する */
	cmp_size = required_dwords;
	if (cmp_size > UINT32_MAX / 2) {
		fputs("size too big!\n", stderr);
		free(calculate_buffer);
		free(results);
		free(result_putterns);
		free(result_all_putterns);
		return 1;
	}
	cmp_size2 = cmp_size * 2;
	cmp_buffer = (uint32_t*)malloc(
		multiply_size(
			multiply_size(sizeof(uint32_t), required_dwords), 5));
	if (cmp_buffer == NULL) {
		perror("failed to allocate compare buffer");
		free(calculate_buffer);
		free(results);
		free(result_putterns);
		free(result_all_putterns);
		return 1;
	}
	cmp_buffer_a = cmp_buffer;
	cmp_buffer_b = cmp_buffer + (size_t)required_dwords * 2;
	cmp_buffer_temp = cmp_buffer + (size_t)required_dwords * 4;

	/* 計算を行う */
	/* サイコロの目の最大値を全探索する */
	for (i = 0; i < dice_max_output; i++) {
		init_dp<<<dice_max_sum / THREADS_ALIGN, THREADS_ALIGN>>>
			(calculate_src, required_dwords);
		cudaDeviceSynchronize();

		/* DP */
		run_dp<<<dice_max_sum / THREADS_ALIGN, THREADS_ALIGN,
			sizeof(uint32_t) * THREADS_ALIGN * required_dwords>>>
			(calculate_src, calculate_dst, required_dwords, i);
		/* サイコロの数繰り返す */
		for (j = 0; j < dice_max_num; j++) {
			/* パターン数を記録する */
			current_result = &results[result_count++];
			current_result->dice_num = j + 1;
			current_result->dice_max_output = i + 1;
			current_result->puttern_count =
				&result_putterns
					[(i * (size_t)dice_max_num + j) * required_dwords];
			current_result->all_puttern_count =
				&result_all_putterns
					[(i * (size_t)dice_max_num + j) * required_dwords];
			if (j == 0) {
				current_result->all_puttern_count[0] = i + 1;
				for (k = 1; k < required_dwords; k++) {
					current_result->all_puttern_count[k] = 0;
				}
			} else {
				multiply(current_result->all_puttern_count,
					previous_all_putterns, i + 1, required_dwords);
			}
			previous_all_putterns = current_result->all_puttern_count;
			cudaDeviceSynchronize();
			cudaMemcpy(current_result->puttern_count,
				&calculate_dst[target_value * (size_t)required_dwords],
				sizeof(*current_result->puttern_count) * required_dwords,
				cudaMemcpyDeviceToHost);
			if (j + 1 < dice_max_num) {
				/* バッファを入れ替える */
				tmp = calculate_src;
				calculate_src = calculate_dst;
				calculate_dst = tmp;
				/* 次の計算を行う */
				run_dp<<<dice_max_sum / THREADS_ALIGN, THREADS_ALIGN,
					sizeof(uint32_t) * THREADS_ALIGN * required_dwords>>>
					(calculate_src, calculate_dst, required_dwords, i);
			}
			/* バッファの挿入ソートを行う */
			for (k = result_count - 1; k > 0; k--) {
				if (result_data_cmp(&results[k - 1], &results[k]) > 0) {
					struct result_data tmp = results[k - 1];
					results[k - 1] = results[k];
					results[k] = tmp;
				} else {
					break;
				}
			}
			/* 表示範囲からあふれた無駄なデータを削る */
			if (result_count > output_num) result_count = output_num;
		}
	}
	cudaFree(calculate_buffer);
	free(cmp_buffer);

	/* 上位を出力する */
	for (i = 0; i < dice_max_sum && i < result_count; i++) {
		printf("%"PRIu32"d%"PRIu32" -> %g (",
			results[i].dice_num, results[i].dice_max_output,
			num_to_double(results[i].puttern_count, required_dwords) /
				num_to_double(results[i].all_puttern_count, required_dwords));
		print_num_and_destroy(results[i].puttern_count, required_dwords);
		printf(" / ");
		print_num_and_destroy(results[i].all_puttern_count, required_dwords);
		printf(")\n");
	}

	free(results);
	free(result_putterns);
	free(result_all_putterns);
	return 0;
}
