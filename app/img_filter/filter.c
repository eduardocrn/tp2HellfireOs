#include <hellfire.h>
#include <noc.h>
#include <math.h>
#include "image.h"

#define MASTER 0
#define BUILDER 1

#define GAUSIAN 0
#define SOBEL 1

#define KEEP_ALIVE 0
#define KILL 3

#define BLOCK_SIZE 36
#define MSG_SIZE 1500

uint8_t gausian(uint8_t buffer[5][5]){
    int32_t sum = 0, mpixel;
    uint8_t i, j;
    int16_t kernel[5][5] =    {    {2, 4, 5, 4, 2},
                    {4, 9, 12, 9, 4},
                    {5, 12, 15, 12, 5},
                    {4, 9, 12, 9, 4},
                    {2, 4, 5, 4, 2}
                };
    for (i = 0; i < 5; i++)
        for (j = 0; j < 5; j++)
            sum += ((int32_t)buffer[i][j] * (int32_t)kernel[i][j]);
    mpixel = (int32_t)(sum / 159);
    return (uint8_t)mpixel;
}

uint32_t isqrt(uint32_t a){
    uint32_t i, rem = 0, root = 0, divisor = 0;
    for (i = 0; i < 16; i++){
        root <<= 1;
        rem = ((rem << 2) + (a >> 30));
        a <<= 2;
        divisor = (root << 1) + 1;
        if (divisor <= rem){
            rem -= divisor;
            root++;
        }
    }
    return root;
}

uint8_t sobel(uint8_t buffer[3][3]){
    int32_t sum = 0, gx = 0, gy = 0;
    uint8_t i, j;
    int16_t kernelx[3][3] =    {    {-1, 0, 1},
                    {-2, 0, 2},
                    {-1, 0, 1},
                };
    int16_t kernely[3][3] =    {    {-1, -2, -1},
                    {0, 0, 0},
                    {1, 2, 1},
                };
    for (i = 0; i < 3; i++){
        for (j = 0; j < 3; j++){
            gx += ((int32_t)buffer[i][j] * (int32_t)kernelx[i][j]);
            gy += ((int32_t)buffer[i][j] * (int32_t)kernely[i][j]);
        }
    }
    sum = isqrt(gy * gy + gx * gx);
    if (sum > 255) sum = 255;
    if (sum < 0) sum = 0;
    return (uint8_t)sum;
}

struct Data{
    // flag, see above
    uint8_t flag;  
    // block order
    uint16_t block_line;
    uint16_t block_column; 
    // pixel block
    uint8_t pixels[BLOCK_SIZE][BLOCK_SIZE]; 
};

union Package{
    struct Data data;
    int8_t raw_data[MSG_SIZE];
};


void init_block(uint8_t block[BLOCK_SIZE][BLOCK_SIZE]){
    uint16_t i, j;
    for (i = 0; i < BLOCK_SIZE; i++){
        for (j = 0; j < BLOCK_SIZE; j++){
            block[i][j] = 255;
        }
    }
}

void do_sobel(uint8_t block[BLOCK_SIZE][BLOCK_SIZE], uint8_t result[BLOCK_SIZE][BLOCK_SIZE]){
    uint16_t i,j,k,l, line, column;
    uint8_t buffer[3][3];

    for (i=0; i < BLOCK_SIZE; i++){
        if (i > 0 && i < BLOCK_SIZE-1){
            for (j=0; j< BLOCK_SIZE; j++){
                if (j > 0 && j < BLOCK_SIZE-1){
                    for (k = 0; k < 3; k++){
                        for(l=0; l<3; l++){
                            line = (i-1) + k; //1 left border + 1 + 1 right border
                            column = (j-1) + l;
                            buffer[k][l] = block[line][column];
                        }
                    }
                    result[i][j] = sobel(buffer);
                }else{
                    result[i][j] = block[i][j];
                }
            }
        }
    }
}

void do_gausian(uint8_t block[BLOCK_SIZE][BLOCK_SIZE], uint8_t result[BLOCK_SIZE][BLOCK_SIZE]){
    uint16_t i,j,k,l, line, column;
    uint8_t buffer[5][5];

    for (i=0; i < BLOCK_SIZE; i++){
        if (i > 1 || i < BLOCK_SIZE-2){
            for (j=0; j< BLOCK_SIZE; j++){
                if (j > 1|| j < BLOCK_SIZE-2){
                    for (k = 0; k < 5; k++){
                        for(l=0; l < 5; l++){
                            line = (i-2) + k; //2 left border + 1 + 2 right border
                            column = (j-2) + l;
                            buffer[k][l] = block[line][column];
                        }
                    }
                    result[i][j] = gausian(buffer);
                }else{
                    result[i][j] = block[i][j];
                }
            }
        }
    }
}

/**
 Populate a block, leaving border of size 2
*/
void create_block(uint8_t matrix[height][width], uint16_t block_line, uint16_t block_column, uint8_t block[BLOCK_SIZE][BLOCK_SIZE]){
    uint16_t start_line = (block_line * (BLOCK_SIZE-4)), line, column;
    uint16_t start_column =  (block_column * (BLOCK_SIZE-4));
    uint8_t i = 2, j;
    //printf("Creating block with start_line %d and start_column %d\n", start_line, start_column);
    while (i < BLOCK_SIZE-2){
        j = 2;
        while (j < BLOCK_SIZE-2){
            line = i+start_line-2;
            column = j+start_column-2;
            if ((line < height) && (column < width)){
                block[i][j] = matrix[line][column];
            }
            j++;
        }
        i++;
    }

}

void master(void){
    int32_t block_lines, block_columns, total_blocks, val, k, l, matrix_pos = 0, next;
    uint8_t max_workes;
    uint16_t i,j, i_line = 0, j_column = 0;

    // int8_t buf[MSG_SIZE];

    // union Package poison_package;

    delay_ms(200);
    
    max_workes = (hf_ncores()-2); // ignoring master and builder

    next = 2;

    block_lines = (height / (BLOCK_SIZE))+2;
    block_columns = (width / (BLOCK_SIZE))+2;

	total_blocks = block_lines * block_columns;

    printf("\n[MASTER] Iniciando... cpu - %d", hf_cpuid());
    printf("\nBlock lines: %d - Block Columns: %d", block_lines-1, block_columns-1);

    if (hf_comm_create(hf_selfid(), 1000, 0))
        panic(0xff);

    // Original image in 2D format
    uint8_t matrix[height][width];

    for (i = 0; i < height; i++){
        for (j = 0; j < width; j++){
            matrix[i][j] = image[matrix_pos];
            matrix_pos++;
        }
    }

    // Block to be processed in the slave cpu
    uint8_t block[BLOCK_SIZE][BLOCK_SIZE];

	while (total_blocks > 0) {

			create_block(matrix, i_line, j_column, block);

            /* Prepare slave message */
            union Package package;
            package.data.flag = GAUSIAN;
            package.data.block_line = i_line;
            package.data.block_column = j_column;
            for (k = 0; k <BLOCK_SIZE; k++)
                for (l = 0; l < BLOCK_SIZE; l++)
                    package.data.pixels[k][l] = block[k][l];

			val = hf_sendack(next, 5000, package.raw_data, MSG_SIZE, next, 700);

            delay_ms(50); //delay 50 or payback timeout 700

            if (val) {
                printf("\nError sending the message to slave. Val = %d \n", val);
			}	

			if (j_column < block_columns) {
				j_column++;
			}
			else {
				i_line++;
				j_column = 0;
            }
            
            total_blocks--;
            next++;
            if (next >= max_workes)
                next = 2;

            delay_ms(20);
	}

	// printf("\nEnviando KILL");
	// // matar escravos
	// poison_package.data.flag = KILL;
    // for (i = 1; i < hf_ncores(); i++) {
	// 	hf_sendack(i, 5000, poison_package.raw_data, MSG_SIZE, i, 500);
	// }

	printf("\n\nend of processing!\n");

    
    hf_kill(hf_selfid());
}

void builder(void) {
    int32_t block_lines, block_columns, total_blocks, val, k, i, j;
    uint16_t cpu, src_port, size;
    int8_t buf[MSG_SIZE];
    int16_t start_line, start_column, line, column, ch;
    union Package package;

    block_lines = (height / (BLOCK_SIZE))+2;
    block_columns = (width / (BLOCK_SIZE))+2;

    total_blocks = block_lines * block_columns;
    
    // Original image in 2D format
    uint8_t matrix[height][width];

    if (hf_comm_create(hf_selfid(), 6000, 0))
        panic(0xff);

    while (total_blocks > 0) {
        // verifica recprob para ver se tem algo a receber
        ch = hf_recvprobe();
        if(ch >= 0 ) {
            // recebe dado do escravo
            val = hf_recvack(&cpu, &src_port, buf, &size, ch);
            if (val){
                printf("\n[slave %d] hf_recvack(): error %d\n", cpu, val);
            } else {
                printf("\n[slave %d] Receive block %d,%d \n", hf_cpuid(), package.data.block_line, package.data.block_column);

                // decode buffer into a union to retrieve struct
                for (i = 0; i < MSG_SIZE; i++)
                    package.raw_data[i] = buf[i];

                // receive pixels, ignoring border and put into main matrix
                start_line = package.data.block_line * (BLOCK_SIZE-4);
                start_column = package.data.block_column * (BLOCK_SIZE-4);                    
                for (i = 2; i < BLOCK_SIZE-2; i++){
                    for (j =2; j < BLOCK_SIZE-2; j++){
                        line = start_line + i - 2;
                        column = start_column + j - 2;
                        if (line < height && column < width)
                            matrix[line][column] = package.data.pixels[i][j];                            
                    }
                }
                
                // atualiza o total_blocks
                total_blocks--;

                printf("\n Remaning blocks %d", total_blocks);
            }

            continue;
        }
    }

    // printf("\nEnviando KILL");
	// // matar escravos
	// poison_package.data.flag = KILL;
    // for (i = 1; i < hf_ncores(); i++) {
	// 	hf_sendack(i, 5000, poison_package.raw_data, MSG_SIZE, i, 500);
	// }


    printf("{CODE}\n");
    printf("\n\nint32_t width = %d, height = %d;\n", width, height);
    printf("uint8_t image[] = {\n");
    k = 0;
    for (i = 0; i < height; i++){
        for (j = 0; j < width; j++){
            printf("0x%x", matrix[i][j]);
            if ((i < height-1) || (j < width-1)) printf(", ");
            if ((++k % 16) == 0) printf("\n");
        }
    }
    printf("};\n");

    printf("{CODE}\n");
}

/**
    Receive pixel vector and applies the requested filter
    Buffer format: [task, position, [pixel_vector]]
    if task is KILL, then this process must die
*/
void slave(void){
    uint8_t keep_alive = 1, result_block[BLOCK_SIZE][BLOCK_SIZE];
    uint16_t cpu, src_port, size, i;
    int16_t val ,ch;
    int8_t buf[MSG_SIZE];
    union Package package;
    init_block(result_block);

    delay_ms(50);
    if (hf_comm_create(hf_selfid(), 5000, 0))
        panic(0xff);

    while (keep_alive == 1){
        ch = hf_recvprobe();
        if (ch >= 0) {
            val = hf_recvack(&cpu, &src_port, buf, &size, ch);
            if (val){
                printf("\n[slave %d] hf_recvack(): error %d\n", cpu, val);
            } else {
                // decode buffer into a union to retrieve struct
                for (i = 0; i < MSG_SIZE; i++){
                    package.raw_data[i] = buf[i];
                }
                printf("\n[slave %d] Receive block %d,%d \n", hf_cpuid(), package.data.block_line, package.data.block_column);
                if (package.data.flag == KILL){
                    keep_alive = 0;
                } else {
                    // package.data.flag = KEEP_ALIVE; 

                    do_gausian(package.data.pixels, result_block);
                    do_sobel(result_block, package.data.pixels);
                    
					delay_ms(random() % 10 + 50);
                    val = hf_sendack(BUILDER, 6000, package.raw_data, MSG_SIZE, hf_cpuid(), 500);

                    if (val)
                        printf("[slave %d] Error sending the message to builder. Val = %d \n", hf_cpuid(), val);
                }
            }
        }
    }

    hf_kill(hf_selfid());
}

void app_main(void) {
    if (hf_cpuid() == MASTER){
        printf("\nCriando Mestre...\n");
        hf_spawn(master, 0, 0, 0, "master", 200000);
    } else if (hf_cpuid() == BUILDER) {
        printf("\nCriando Builder...\n");
        hf_spawn(builder, 0, 0, 0, "builder", 200000);
    } else{
        printf("\nCriando escravo %d...\n", hf_cpuid());
        hf_spawn(slave, 0, 0, 0, "slave", 100000);
    }
}
