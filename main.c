#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#define MAX 8
#define HI(num)	(((num) & 0x0000FF00) << 8)
#define LO(num)	((num) & 0x000000FF)

void printer(double matrix[MAX][MAX]){
    int i,j;
    printf("\n");
    for(i=0;i<MAX;i++){
        for(j=0;j<MAX;j++){
            printf("%f ", matrix[i][j] );
        }
        printf("\n");
    }
    printf("\n");
}

int block_creation(double block[MAX][MAX], double qtrd[MAX][MAX]){
    int n,k,l;
    int zeros=0;
    //decrase 128 from each value

    for(k=0;k<MAX;k++){
      for(n=0;n<MAX;n++) block[k][n] = block[k][n]-128;
    }
    //printf("Leveled off Block:\n");
    //printer(block);

    //matrix DCT
    double dct[MAX][MAX],transformed[MAX][MAX],tdct[MAX][MAX];

    for(n=0;n<MAX;n++) dct[0][n]=1/sqrt(MAX);

    for(k=1;k<MAX;k++)
      for(n=0;n<MAX;n++) dct[k][n] = sqrt(2.0/MAX) * cos(k*M_PI/(2.0*MAX)*(2.0*n+1));
    //printf("DCT Matrix:\n");
    //printer(dct);

    //transposta de DCT
    for (n = 0; n < MAX; n++) {
        for (k =0; k < MAX; k++) {
            tdct[k][n]=dct[n][k];
        }
    }
    //printer(tdct);

    // DCT * Block

    for(n=0;n<MAX;n++){
        for(k=0;k<MAX;k++){
            transformed[n][k]=0;
            for(l=0;l<MAX;l++) transformed[n][k] += dct[n][l] * block[l][k] ;
        }
    }

    // Block * DCT(transposta)
    for(n=0;n<MAX;n++){
        for(k=0;k<MAX;k++){
            block[n][k]=0;
            for(l=0;l<MAX;l++) block[n][k] += transformed[n][l] * tdct[l][k];
        }
    }
    //printf("DCT x Block x DCT(transposed):\n");
    //printer(block);

    //decompressor block_decompressor(block);

    //ajusting the desired quality percentage


    for(n=0;n<MAX;n++){
        for(k=0;k<MAX;k++){
            block[n][k] = (int)round( block[n][k] / qtrd[n][k]);
        }
    }

    //zigzaged stream

    //high part
    int zigzag[64];
    int nx=0,dir=1;
    l=0;
    do{
        if(dir==1){
            k=0;
            for(n=nx;n>=0;n--){
                zigzag[l]=block[n][k];
                if (k<nx){
                    k++;
                }
                l++;
            }
            dir=0;
        }else{
            k=nx;
            for(n=0;n<=nx;n++){
                zigzag[l]=block[n][k];
                if (k>0){
                    k--;
                }
                l++;
            }
            dir=1;
        }
        nx++;
    }while(nx<MAX);
    //low part

    l=63;dir=0;nx=7;
    do{
        if(dir==1){
            k=nx;
            for(n=MAX-1;n>=nx;n--){
                zigzag[l]=block[n][k];
                if (k<MAX){
                    k++;
                }
                l--;
            }
            dir=0;
        }else{
            k=MAX-1;
            for(n=nx;n<=MAX;n++){
                zigzag[l]=block[n][k];
                if (k>0){
                    k--;
                }
                l--;
            }
            dir=1;
        }
        nx--;
    }while(nx>0);

    n=63;
    while((zigzag[n]==0) && (n>=0) ){
        zeros++;
        n--;
    }

    return zeros;
}

void block_decompressor(double block[MAX][MAX], double qtrd[MAX][MAX]){
    double transformed[MAX][MAX];
    int n,k,l;
    double dct[MAX][MAX],tdct[MAX][MAX];

    for(n=0;n<MAX;n++) dct[0][n]=1/sqrt(MAX);

    for(k=1;k<MAX;k++)
      for(n=0;n<MAX;n++) dct[k][n] = sqrt(2.0/MAX) * cos(k*M_PI/(2.0*MAX)*(2.0*n+1));


    //transposta de DCT
    for (n = 0; n < MAX; n++) {
        for (k =0; k < MAX; k++) {
            tdct[k][n]=dct[n][k];
        }
    }

    //R(i,j) = Quantization(n,k) * Block(n*k)


    for(n=0;n<MAX;n++){
        for(k=0;k<MAX;k++){
            block[n][k] = block[n][k] * qtrd[n][k];
        }
    }

    // DCT(transposta) * Block

    for(n=0;n<MAX;n++){
        for(k=0;k<MAX;k++){
            transformed[n][k]=0;
            for(l=0;l<MAX;l++) transformed[n][k] += tdct[n][l] * block[l][k] ;
        }
    }

    // Block * DCT
    for(n=0;n<MAX;n++){
        for(k=0;k<MAX;k++){
            block[n][k]=0;
            for(l=0;l<MAX;l++) block[n][k] += transformed[n][l] * dct[l][k];
        }
    }

    //round +128

    for(n=0;n<MAX;n++){
        for(k=0;k<MAX;k++){
            block[n][k] = round(block[n][k])+128;
        }
    }
}

typedef struct PGMatrix {
    int row;
    int col;
    int max_gray;
    int **matrix;
    char Pn[3];
} PGMData;

int **allocate_dynamic_matrix(int row, int col){
    int **ret_val;
    int i;
    ret_val = (int **)malloc(sizeof(int *) * row);
    if (ret_val == NULL) {
        perror("memory allocation failure");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < row; ++i) {
        ret_val[i] = (int *)malloc(sizeof(int) * col);
        if (ret_val[i] == NULL) {
            perror("memory allocation failure");
            exit(EXIT_FAILURE);
        }
    }
    return ret_val;
}

void deallocate_dynamic_matrix(int **matrix, int row){
    int i;
    for (i = 0; i < row; ++i) {
        free(matrix[i]);
    }
    free(matrix);
}
void SkipComments(FILE *fp){
    int ch;
    char line[100];
    while ((ch = fgetc(fp)) != EOF && isspace(ch)) {
        ;
    }
    if (ch == '#') {
        fgets(line, sizeof(line), fp);
        SkipComments(fp);
    } else {
        fseek(fp, -1, SEEK_CUR);
    }
}

/*for reading:*/
PGMData* readPGM(const char *file_name, PGMData *data){
    FILE *pgmFile;
    char version[3];
    int i, j;
    int lo, hi;
    pgmFile = fopen(file_name, "rb");
    if (pgmFile == NULL) {
        perror("cannot open file to read");
        exit(EXIT_FAILURE);
    }
    fgets(version, sizeof(version), pgmFile);
    strcpy(data->Pn,version);
    if (strcmp(version, "P5")) {
        fprintf(stderr, "Wrong file type!\n");
        exit(EXIT_FAILURE);
    }
    SkipComments(pgmFile);
    fscanf(pgmFile, "%d", &data->col);
    SkipComments(pgmFile);
    fscanf(pgmFile, "%d", &data->row);
    SkipComments(pgmFile);
    fscanf(pgmFile, "%d", &data->max_gray);
    fgetc(pgmFile);

    data->matrix = allocate_dynamic_matrix(data->row, data->col);
    if (data->max_gray > 255) {
        for (i = 0; i < data->row; ++i) {
            for (j = 0; j < data->col; ++j) {
                hi = fgetc(pgmFile);
                lo = fgetc(pgmFile);
                data->matrix[i][j] = (hi << 8) + lo;
            }
        }
    }
    else {
        for (i = 0; i < data->row; ++i) {
            for (j = 0; j < data->col; ++j) {
                lo = fgetc(pgmFile);
                data->matrix[i][j] = lo;
            }
        }
    }

    fclose(pgmFile);
    return data;

}
void initdata(PGMData * data){
    data->row=0;
    data->col=0;
    data->max_gray=0;
    data->Pn[0]='P';
    data->Pn[1]='2';
}

void writePGM(const char *filename, const PGMData *data)
{
    FILE *pgmFile;
    int i, j;
    int hi, lo;

    pgmFile = fopen(filename, "wb");
    if (pgmFile == NULL) {
        perror("cannot open file to write");
        exit(EXIT_FAILURE);
    }
    fprintf(pgmFile, "P5 ");
    fprintf(pgmFile, "%d %d ", data->col, data->row);
    fprintf(pgmFile, "%d ", data->max_gray);

    if (data->max_gray > 255) {
        for (i = 0; i < data->row; ++i) {
            for (j = 0; j < data->col; ++j) {
                hi = HI(data->matrix[i][j]);
                lo = LO(data->matrix[i][j]);
                fputc(hi, pgmFile);
                fputc(lo, pgmFile);
            }

        }
    }
    else {
        for (i = 0; i < data->row; ++i) {
            for (j = 0; j < data->col; ++j) {
                lo = LO(data->matrix[i][j]);
                fputc(lo, pgmFile);
            }
        }
    }

    fclose(pgmFile);
    deallocate_dynamic_matrix(data->matrix, data->row);
}

int main(int argc, const char * argv[])
{
    int zeros=0;
    PGMData data,data2;
    initdata(&data);initdata(&data2);
    readPGM("lena.pgm",&data);
    readPGM("lena.pgm",&data2);
    int i=0, j=0,ik,jk,temp=0;
    double matrix[MAX][MAX];

    double qt50[MAX][MAX]={
        {16, 11, 10, 16, 24, 40, 51, 61},
        {12, 12, 14, 19, 26, 58, 60, 55},
        {14, 13, 16, 24, 40, 57, 69, 56},
        {14, 17, 22, 29, 51, 87, 80, 62},
        {18, 22, 37, 56, 68, 109, 103, 77},
        {24, 35, 55, 64, 81, 104, 113, 92},
        {49, 64, 78, 87, 103, 121, 120, 101},
        {72, 92, 95, 98, 112, 100, 103, 99}
    };

    double qt90[MAX][MAX]={
        {80, 60, 50, 80, 120, 200, 255, 255},
        {55, 60, 70, 95, 130, 255, 255, 255},
        {70, 65, 80, 120, 200, 255, 255, 255},
        {70, 85, 110, 145, 255, 255, 255, 255},
        {90, 110, 185, 255, 255, 255, 255, 255},
        {120, 175, 255, 255, 255, 255, 255, 255},
        {245, 255, 255, 255, 255, 255, 255, 255},
        {255, 255, 255, 255, 255, 255, 255, 255}
    };

    double qtrd[MAX][MAX];

    for(i=0;i<MAX;i++){
        for(j=0;j<MAX;j++){
            qtrd[i][j]=255;
        }
    }
    qtrd[0][0]=80;


    for(i=0;i<data.col;i+=8){
        for(j=0;j<data.row;j+=8){

            for(ik=0;ik<MAX;ik++){
                for(jk=0;jk<MAX;jk++){
                    matrix[ik][jk]=data.matrix[ik+i][jk+j];
                }
            }
            zeros+=block_creation(&matrix,qtrd);
            block_decompressor(&matrix,qtrd);

            for(ik=0;ik<MAX;ik++){
                for(jk=0;jk<MAX;jk++){
                    data2.matrix[ik+i][jk+j]=(int)matrix[ik][jk];
                }
            }

        }
    }


    printf("\nVALUE\n%d\n",zeros);
    writePGM("l.pgm",&data2);
    return 0;
}

