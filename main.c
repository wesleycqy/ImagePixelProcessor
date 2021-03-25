// Wesley Chong CIS 416
// gcc -O1 main.c -o pixels; ./pixels

// TO SEE MAXIMUM EFFICIENCY IMPROVEMENTS (UP TO > 6x) FOR ROTATION, PLEASE RUN IN OPTIMIZATION LEVEL 1 (-O1)

// TO SEE MAXIMUM EFFICIENCY IMPROVEMENTS (UP TO > 2x) FOR SMOOTH, PLEASE RUN IN OPTIMIZATION LEVEL 0 (-O0)

#define N (1<<10) 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

typedef unsigned short pixel;

void pack_pixel(float r, float g, float b, pixel *p)
{
    *p = (((int) (r * 31) & 0x1F) << 11) +
         (((int) (g * 63) & 0x3F) << 5) +
          ((int) (b * 31) & 0x1F);
}
void unpack_pixel(pixel p, float *r, float *g, float *b) 
{
  *r = (p>>11) / 31.0f;
  *g = ((p>>5) & 0b111111) / 63.0f;
  *b = (p & 0x1F) / 31.0f;
}

void smooth_pack_pixel(float r, float g, float b, pixel *p) // uses int instead of float
{
    *p = ((int) r  << 11) + (((int) g  & 0x3F) << 5) + ((int) b  & 0x1F);
}

void smooth_unpack_pixel(pixel p, int *r, int *g, int *b) // uses int instead of float
{
  *r = ( (p >> 11) ); 
  *g = ( ((p >> 5) & 0x3F) );
  *b = ( (p & 0x1F) );
}

void pixel_to_gray(pixel p, pixel *q)
{
    float r, g, b, avg;
    unpack_pixel(p, &r, &g, &b);
    avg = (r + g + b) / 3;
    pack_pixel(avg, avg, avg, q);
}

void array_to_gray(pixel *p, pixel *q)
{
    int i, j;
    
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            pixel_to_gray(p[i * N + j], &q[i * N + j]);
}

void pixel_to_sepia(pixel p, pixel *q)
{
    float r, g, b;
    unpack_pixel(p, &r, &g, &b);


    float new_r = r * .393 + g * .769 + b * .189;
    if (new_r >32){new_r = 32;}

    float new_g = r * .349 + g * .686 + b * .168;
    if (new_g >32){new_g = 32;}

    float new_b = r * .272 + g * .534 + b * .131;
    if (new_b >32){new_b = 32;}

    pack_pixel(new_r, new_g, new_b, q);
}

void array_to_sepia(pixel *p, pixel *q)
{
  int i, j;
    
  for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
          pixel_to_sepia(p[i * N + j], &q[i * N + j]);
}

void transpose(pixel *p, pixel *q)    // original
{
    int i, j;
    
    for (i = 0; i < N; i++)
        for (j = 0; j <= i; j++)
        {
            q[i * N + j] = p[j * N + i];
            
            if (j != i)
                q[j * N + i] = p[i * N + j];
        }
}

void transpose2(pixel *p, pixel *q)   // speeds up multiplication for i loop
{
  int i, j;
  
  for (i = 0; i < N; i++){
    int iN = i * N;
    for (j = 0; j <= i; j++)
      {
        q[iN + j] = p[j * N + i];
        if (j != i)
            q[j * N + i] = p[iN + j];
      }
  }
}

void transpose3(pixel *p, pixel *q)   // variable declaration
{
  register unsigned int i, j;
  
  for (i = 0; i < N; i++){
    int iN = i * N;
    for (j = 0; j <= i; j++)
    {
      int jN = j * N;
      q[iN + j] = p[jN + i];
      if (j != i){q[jN + i] = p[iN + j];}
    }
  }
}

void transpose4(pixel *p)             // remove need for second argument
{
  register unsigned int i, j;
  int iN, jN, t;
  
  for (i = 0, iN = 0; i < N; i++, iN+=N){
    for (j = 0, jN = 0; j <= i; j++, jN += N)
    {
      t = p[iN + j];
      p[iN + j] = p[jN + i];
      if (j != i){p[jN + i] = t;}
    }
  }
}


void flip_vertical(pixel *p, pixel *q)      // original
{
    int i, j;
    
    for (i = 0; i <= N / 2; i++)
        for (j = 0; j < N; j++)
        {
            q[(N - i - 1) * N + j] = p[i * N + j];
            
            if (i < N / 2)
                q[i * N + j] = p[(N - i - 1) * N + j];
        }
}

void flip_vertical2(pixel *p, pixel *q)     // speeds up i loop multiplication
{
    int i, j;
    
    for (i = 0; i <= N / 2; i++)
    {
      int iN = i * N;
      int calculation = (N - i - 1) * N;
      for (j = 0; j < N; j++)
      {
          q[calculation + j] = p[iN + j];
          
          if (i < N / 2)
              q[iN + j] = p[calculation + j];
      }
    }
}

void flip_vertical3(pixel *p, pixel *q)     // loop unrolling, variable declaration, & other improvements
{
    register unsigned int i, j;
    
    for (i = 0; i <= N >> 1; i++)
    {
      int iN = i * N;
      int calculation = (N - i - 1) * N;

      for (j = 0; j < N - 1; j+= 2)
      {
          q[calculation + j] = p[iN + j];
          
          if (i < N >> 1){q[iN + j] = p[calculation + j];}

          q[calculation + j + 1] = p[iN + j + 1];
          
          if (i < N >> 1){q[iN + j + 1] = p[calculation + j + 1];}
      }
      if(j % 2 == 1){
        q[calculation + j - 1] = p[iN + j - 1];  
        if (i < N >> 1){q[iN + j - 1] = p[calculation + j - 1];}
      }
    }
}

void flip_horizontal(pixel *p, pixel *q)
{
  for (int j = 0; j < N ; j++) {
    for (int i = j * 5; i < N * (j + 1); i++) {
      q[i] = p[i - 2 * (i % N) + N - 1];
    }
  }
}


void rotate1(pixel *p, pixel *q){
    pixel *t = (pixel *) calloc(N * N, sizeof(pixel));
    transpose(p, t);
    flip_vertical(t, q);
    free(t);
} // original

void rotate2(pixel *p, pixel *q){
  pixel * t = (pixel *) calloc(N * N, sizeof(pixel));
  transpose2(p,t);
  flip_vertical2(t, q);
  free(t);
} // New transpose & New flip_vertical

void rotate3(pixel *p, pixel *q){
  pixel * t = (pixel *) calloc(N * N, sizeof(pixel));
  transpose3(p,t);
  flip_vertical3(t, q);
  free(t);
} // New transpose & New flip_vertical

void rotate4(pixel *p, pixel *q){
  pixel * t = (pixel *) calloc(N * N, sizeof(pixel));
  memcpy(q, p, N * N * sizeof(pixel));
  transpose4(q);
  flip_vertical3(t, q);
} // New transpose only

void average1(pixel *p, int i, int j, pixel *q, int pixelNumber)
{
  // -N-1, -N, -N+1,       
  // -1, 0, +1,
  // N-1, N, N+1
  float r, g, b;
  float count=9;
  float total_r = 0 ;
  float total_g = 0 ;
  float total_b = 0 ;

  if(pixelNumber==0){
    count=4;

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N+1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if((N-1) > pixelNumber >= 1){
    count=6;

    unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N+1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if(pixelNumber==(N - 1)){
    count = 4;

    unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if(pixelNumber==(N * N) - (N)){
    count = 4; 

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
    
    unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if(pixelNumber % N ==0){
    count = 6;
    
    unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;


  }
  else if((pixelNumber + 1) % N ==0){
    count = 6;
    
    unpack_pixel(p[-N - 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N - 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if((N * N - 2) >= pixelNumber && pixelNumber > ((N * N) - N) ){
    count = 6;
    
    unpack_pixel(p[-N - 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if(pixelNumber== (N * N - 1)){
    count = 4;

    unpack_pixel(p[-N - 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else{
    count = 9;

    unpack_pixel(p[-N-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N - 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  
  float new_r = total_r / count;
  float new_g = total_g / count;
  float new_b = total_b / count; 

  pack_pixel(new_r, new_g, new_b, q);
} // original

void average2(pixel *p, int i, int j, pixel *q, int pixelNumber) // prioritize multiplication
{
  // -N-1, -N, -N+1,       
  // -1, 0, +1,
  // N-1, N, N+1
  float r, g, b;
  float count=9;
  float total_r = 0 ;
  float total_g = 0 ;
  float total_b = 0 ;
  int NN = N * N;

  if(pixelNumber==0){
    count=4;

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N+1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if((N-1) > pixelNumber >= 1){
    count=6;

    unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N+1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if(pixelNumber==(N - 1)){
    count = 4;

    unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if(pixelNumber==(N * N) - (N)){
    count = 4; 

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
    
    unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if(pixelNumber % N ==0){
    count = 6;
    
    unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;


  }
  else if((pixelNumber + 1) % N ==0){
    count = 6;
    
    unpack_pixel(p[-N - 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N - 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if((NN - 2) >= pixelNumber && pixelNumber > ((NN) - N) ){
    count = 6;
    
    unpack_pixel(p[-N - 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if(pixelNumber== (NN - 1)){
    count = 4;

    unpack_pixel(p[-N - 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else{
    count = 9;

    unpack_pixel(p[-N-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N - 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    unpack_pixel(p[N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  
  float new_r = total_r / count;
  float new_g = total_g / count;
  float new_b = total_b / count; 

  pack_pixel(new_r, new_g, new_b, q);
}

void average3(pixel *p, int i, int j, pixel *q, int pixelNumber) // use ints when unpacking & packing pixels
{
  // -N-1, -N, -N+1,       
  // -1, 0, +1,
  // N-1, N, N+1
  int r, g, b;
  float count=9;
  float total_r = 0 ;
  float total_g = 0 ;
  float total_b = 0 ;
  int NN = N * N;

  if(pixelNumber==0){
    count=4;

    smooth_unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[N+1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if((N-1) > pixelNumber >= 1){
    count=6;

    smooth_unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[N-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[N+1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if(pixelNumber==(N - 1)){
    count = 4;

    smooth_unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[N-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if(pixelNumber==(N * N) - (N)){
    count = 4; 

    smooth_unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
    
    smooth_unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[-N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if(pixelNumber % N ==0){
    count = 6;
    
    smooth_unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[-N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if((pixelNumber + 1) % N ==0){
    count = 6;
    
    smooth_unpack_pixel(p[-N - 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[N - 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if((NN - 2) >= pixelNumber && pixelNumber > ((NN) - N) ){
    count = 6;
    
    smooth_unpack_pixel(p[-N - 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[-N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else if(pixelNumber== (NN - 1)){
    count = 4;

    smooth_unpack_pixel(p[-N - 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  else{
    count = 9;

    smooth_unpack_pixel(p[-N-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[-N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[-N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[-1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[0], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[N - 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[N], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;

    smooth_unpack_pixel(p[N + 1], &r, &g, &b);
    total_r += r ;
    total_g += g ;
    total_b += b ;
  }
  
  float new_r = total_r / count;
  float new_g = total_g / count;
  float new_b = total_b / count; 

  smooth_pack_pixel(new_r, new_g, new_b, q);
}


void smooth1(pixel *p, pixel *q)
{
    int i, j;
    float r, g, b;
    int count = 0;

    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++)
        {
          average1(&p[count], i, j, &q[i * N + j], count);
          count++;
        }
    }    
} // original

void smooth2(pixel *p, pixel *q) 
{
    register unsigned int i, j;
    float r, g, b;
    register unsigned int count = 0;
    int iN = i * N;


    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++)
        {
          average2(&p[count], i, j, &q[iN + j], count);
          count++;
        }
    }    
} // change loop counter to registers & speed up multiplication

void smooth3(pixel *p, pixel *q) // use ints while packing/unpacking pixels in avg function
{
    register unsigned int i, j;
    float r, g, b;
    register unsigned int count = 0;


    for (i = 0; i < N; i++){
        for (j = 0; j < N; j++)
        {
          average3(&p[count], i, j, &q[i * N + j], count);
          count++;
        }
    }    
}

void print_array(pixel *p)
{
    int i, j;
    float r, g, b;
    
    if (N > 6)
    {
        //printf("TOO BIG TO PRINT!\n");
        return;
    }
    
    printf("\n");

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            unpack_pixel(*(p + i * N + j), &r, &g, &b);

            if (r == 0 && g == 0 && b == 0)
                printf("-----------  ");
            else
                printf("%.1f,%.1f,%.1f  ", r, g, b);
        }
            
        printf("\n");
    }
}

int main(int argc, char **argv)
{
    long start;
    pixel *p = (pixel *) calloc(N * N, sizeof(pixel));
    pixel *q = (pixel *) calloc(N * N, sizeof(pixel));
   

    pack_pixel(0.1, 0.1, 0.1, p);
    pack_pixel(0.2, 0.2, 0.2, &p[1]);
    pack_pixel(0.3, 0.3, 0.3, &p[2]);
    pack_pixel(0.4, 0.4, 0.4, &p[N]);
    pack_pixel(0.9, 0.9, 0.9, &p[N * N - 1]);
   
   // TO SEE MAXIMUM EFFICIENCY IMPROVEMENTS (UP TO > 6x) FOR ROTATION, PLEASE RUN IN OPTIMIZATION LEVEL 1 (-O1)

    start = clock();
    rotate1(p, q);
    printf("\nRotates:\nrotate1: %.3f seconds\n",(clock() - start) / (double) CLOCKS_PER_SEC);

    start = clock();
    rotate2(p, q);
    printf("rotate2: %.3f seconds\n",(clock() - start) / (double) CLOCKS_PER_SEC);

    start = clock();
    rotate3(p, q);
    printf("rotate3: %.3f seconds\n",(clock() - start) / (double) CLOCKS_PER_SEC);

    start = clock();
    rotate4(p, q);
    printf("rotate4: %.3f seconds\n\n",(clock() - start) / (double) CLOCKS_PER_SEC);

  // TO SEE MAXIMUM EFFICIENCY IMPROVEMENTS (UP TO > 2x) FOR SMOOTH, PLEASE RUN IN OPTIMIZATION LEVEL 0 (-O0)

    // start = clock();
    // smooth1(p, q);
    // printf("smooth1: %.3f seconds\n\n",(clock() - start) / (double) CLOCKS_PER_SEC);

    // start = clock();
    // smooth2(p, q);
    // printf("smooth2: %.3f seconds\n\n",(clock() - start) / (double) CLOCKS_PER_SEC);

    // start = clock();
    // smooth3(p, q);
    // printf("smooth3: %.3f seconds\n\n",(clock() - start) / (double) CLOCKS_PER_SEC);

    print_array(p);
    print_array(q);
    
    return 0;
}