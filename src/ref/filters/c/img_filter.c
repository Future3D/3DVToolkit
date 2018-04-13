/**
 * @file    img_filter.c
 * @author  YRD
 * @version 1.0
 * @date    2016-9-20
 * @brief   图像滤波运算函数
 * @details 本文件包括了时域和空间域的滤波算法
*/


#ifdef WIN32
#pragma warning (disable:4996)
#endif

#include "../comm/api.h"

#include "img_const.h"
#include "img_api.h"
#include "img_algo.h"
#include "img_filter.h"
#include <string.h>
#include <math.h>
#include <stdint.h>

/** 
 * @fn              float *img_fir_cross(float *img_out, float *img_in, float *coff)
 * @details         使用十字滤波模板的图像滤波，滤波器模板如下
 *                             0(p)
 *                   1(p+W-1)  2(p+W)  3(p+W+1)
 *                             4(p+2W)
 *                  注意：滤波输出图像的最外圈边沿（1层像素）是无效数据
 * @param [in]      float *img_in：指针，指向待滤波图像
 * @param [in]      float *coff：指针，指向5个滤波系数
 * @param [out]     float *img_out：指针，指向的空间存放图像运算结果
 * @retval          float *：和img_out相同
 */ 
float *img_fir_cross(float *img_out, float *img_in, float *coff)
{
    float s;
    
    float *p0=img_in+1;
    float *p1=p0+IMG_WID-1;
    float *p2=p0+IMG_WID  ;
    float *p3=p0+IMG_WID+1;
    float *p4=p0+IMG_WID+IMG_WID;
    
    float *q=img_out+IMG_WID+1,*q_end=img_out+IMG_WID*(IMG_HGT-1)-1;
    
    float c0=*(coff+0), c1=*(coff+1), c2=*(coff+2), c3=*(coff+3), c4=*(coff+4);
    
    for (;q<q_end;q++,p0++,p1++,p2++,p3++,p4++)
    {
        s =(*p0)*c0;
        s+=(*p1)*c1;
        s+=(*p2)*c2;
        s+=(*p3)*c3;
        s+=(*p4)*c4; 

        *q=s;
    }
    return img_out;
}


/** 
 * @fn              float *img_fir_cross_sa(float *img_out, float *img_in, float *coff)
 * @details         使用十字滤波模板的图像滤波（原址操作），功能同img_fir_cross，但原址运算实现
 *                  注意：滤波输出图像的最外圈边沿（1层像素）是无效数据
 * @param [inout]   float *img_inout：指针，指向待滤波图和图像运算结果
 * @param [in]      float *coff：指针，指向5个滤波系数
 * @param [inout]   环形缓冲器，存放3行数据
 * @retval          float *：和img_inout相同
 */ 
float *img_fir_cross_sa(float *img_inout, float *coff, struct ring_buf_f32_s *rbuf)
{
    float s;
    
    float *p0=img_inout+1;
    float *p1=p0+IMG_WID-1;
    float *p2=p0+IMG_WID  ;
    float *p3=p0+IMG_WID+1;
    float *p4=p0+IMG_WID+IMG_WID;
    
    float *q =img_inout+IMG_WID+1-3*IMG_WID,*q_end=img_inout+IMG_WID*(IMG_HGT-1)-1-3*IMG_WID;
    
    float c0=*(coff+0), c1=*(coff+1), c2=*(coff+2), c3=*(coff+3), c4=*(coff+4);
    
    int n=3*IMG_WID;
    for (;q<q_end;q++,p0++,p1++,p2++,p3++,p4++)
    {
        s =(*p0)*c0;
        s+=(*p1)*c1;
        s+=(*p2)*c2;
        s+=(*p3)*c3;
        s+=(*p4)*c4; 
        
        if (n)
        {
            n--;
            ring_buf_f32_io(rbuf,s);
        }
        else
            *q=ring_buf_f32_io(rbuf,s);
    }

    for (n=0;n<3*IMG_WID;n++,q++)
        *q=ring_buf_f32_io(rbuf,0);

    return img_inout;
}


/** 
 * @fn              float *img_fir_sqr3(float *img_out, float *img_in, float *coff)
 * @details         使用3x3滤波模板的图像滤波，滤波器模板如下
 *                  0(p)    1(p+1)    2(p+2)
 *                  3(p+W)  4(p+W+1)  5(p+W+2)
 *                  6(P+2W) 7(p+2W+1) 8(p+2W+2)
 *                  注意：滤波输出图像的最外圈边沿（1层像素）是无效数据
 * @param [in]      float *img_in：指针，指向待滤波图像
 * @param [in]      float *coff：指针，指向9个滤波系数
 * @param [out]     float *img_out：指针，指向的空间存放图像运算结果
 * @retval          float *：和img_out相同
 */ 
float *img_fir_sqr3(float *img_out, float *img_in, float *coff)
{
    float s;
    float *p0=img_in          , *p1=img_in          +1, *p2=img_in          +2;
    float *p3=img_in+  IMG_WID, *p4=img_in+  IMG_WID+1, *p5=img_in+  IMG_WID+2;
    float *p6=img_in+2*IMG_WID, *p7=img_in+2*IMG_WID+1, *p8=img_in+2*IMG_WID+2;

    float *q=img_out+IMG_WID+1,*q_end=img_out+IMG_WID*(IMG_HGT-1)-1;
    
    float c0=*(coff  ), c1=*(coff+1), c2=*(coff+2);
    float c3=*(coff+3), c4=*(coff+4), c5=*(coff+5);
    float c6=*(coff+6), c7=*(coff+7), c8=*(coff+8);

    for (;q<q_end;q++,p0++,p1++,p2++,p3++,p4++,p5++,p6++,p7++,p8++)
    {
        s =(*p0) * c0;
        s+=(*p1) * c1;
        s+=(*p2) * c2;
        s+=(*p3) * c3;
        s+=(*p4) * c4;
        s+=(*p5) * c5;
        s+=(*p6) * c6;
        s+=(*p7) * c7;
        s+=(*p8) * c8;

        *q=s;
    }

    return img_out;
}


/** 
 * @fn              float *img_fir_sqr3_sa(float *img_inout, float *coff, struct ring_buf_f32_s *rbuf)
 * @details         使用3x3滤波模板的图像滤波（原址操作），功能同img_fir_sqr3，但通过原址操作实现
 *                  注意：滤波输出图像的最外圈边沿（1层像素）是无效数据
 * @param [inout]   float *img_inout：指针，指向待滤波图和图像运算结果
 * @param [in]      float *coff：指针，指向9个滤波系数
 * @param [inout]   ring_buf_f32_s *rbuf环形缓冲器，存放3行数据
 * @retval          float*：和img_inout相同
 */ 
float *img_fir_sqr3_sa(float *img_inout, float *coff, struct ring_buf_f32_s *rbuf)
{
    float s;
    float *p0=img_inout          , *p1=img_inout          +1, *p2=img_inout          +2;
    float *p3=img_inout+  IMG_WID, *p4=img_inout+  IMG_WID+1, *p5=img_inout+  IMG_WID+2;
    float *p6=img_inout+2*IMG_WID, *p7=img_inout+2*IMG_WID+1, *p8=img_inout+2*IMG_WID+2;

    float *q=img_inout+IMG_WID+1-3*IMG_WID,*q_end=img_inout+IMG_WID*(IMG_HGT-1)-1-3*IMG_WID;
    
    int n=3*IMG_WID;

    float c0=*(coff  ), c1=*(coff+1), c2=*(coff+2);
    float c3=*(coff+3), c4=*(coff+4), c5=*(coff+5);
    float c6=*(coff+6), c7=*(coff+7), c8=*(coff+8);

    for (;q<q_end;q++,p0++,p1++,p2++,p3++,p4++,p5++,p6++,p7++,p8++)
    {
        s =(*p0) * c0;
        s+=(*p1) * c1;
        s+=(*p2) * c2;
        s+=(*p3) * c3;
        s+=(*p4) * c4;
        s+=(*p5) * c5;
        s+=(*p6) * c6;
        s+=(*p7) * c7;
        s+=(*p8) * c8;

        if (n)
        {
            n--;
            ring_buf_f32_io(rbuf,s);
        }
        else
            *q=ring_buf_f32_io(rbuf,s);
    }

    for (n=0;n<3*IMG_WID;n++,q++)
        *q=ring_buf_f32_io(rbuf,0);

    return img_inout;
}


/** 
 * @fn              float *img_iir_t(float *img_out,float *img_in, float alpha)
 * @details         1阶IIR图像序列的时间滤波，使用有损积分器结构
 * @param [in]      float *img_in：指针，指向待滤波图像
 * @param [in]      float alpha：滤波系数（遗忘因子）0~1，越接近1，滤波器带宽越小
 * @param [inout]   float *img_inout：指针，指向空间存放先前滤波结果和新的滤波结果
 * @retval          float *：和img_inout相同
 */
float *img_iir_t(float *img_inout,float *img_in, float alpha)
{
    float *p=img_in, *p_end=img_in+IMG_WID*IMG_HGT;
    float *q=img_inout;
    for (;p<p_end;p++,q++)
        (*q)=(*q)*(float)alpha+(float)(1.0-alpha)*(*p);
    return img_inout;
}


/** 
 * @fn              float *img_fir3_t_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2, float *coff)
 * @details         图像序列的FIR时间滤波，使用前2帧和当前帧数据
 * @param [in]      float *img0，img1，img2为历史图像帧(指针)，img0对应最老图像，img2对应最新图像
 * @param [in]      float *coff：指针，指向滤波加权系数数组，*coff对应img_in0，*(coff+1)对应img_in1,*(coff+2)对应img_in2
 * @param [out]     float *img_out：指针，指向空间存放滤波结果
 * @retval          float *：和img_out相同
 */
float *img_fir3_t_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2, float *coff)
{
    float *p0=img_in0, *p1=img_in1, *p2=img_in2;
    float *q=img_out,*q_end=img_out+IMG_WID*IMG_HGT;
    float c0=*(coff),c1=*(coff+1),c2=*(coff+2);

    float s;
    for (;q<q_end;p0++,p1++,p2++,q++)
    {
        s =(*p0)*c0;
        s+=(*p1)*c1;
        s+=(*p2)*c2;
        *q=s;
    }
    return img_out;
}


/** 
 * @fn              float *img_fir3_t(float *img_out, float *img_in0, float *img_in1, float *img_in2, float *coff)
 * @details         图像序列的FIR时间滤波，使用前2帧和当前帧数据
 * @param [in]      float *img_buf：为历史图像帧(指针)，指向区域连续存放最近2帧图像
 * @param [in]      float *img_in：指针，指向最新输入图像
 * @param [in]      float *coff：指针，指向滤波加权系数数组，*coff对应img_in0，*(coff+1)对应img_in1,*(coff+2)对应img_in2
 * @param [out]     float *img_out：指针，指向的空间存放滤波结果
 * @param [inout]   int *state：指针，指向滤波状态变量，初始值需设为0，指向的内容在运行后被修改
 * @retval          float *：和img_out相同
 */
float *img_fir3_t(float *img_out, float *img_buf, float *img_in, float *coff, int *state)
{
    float *img_buf1=img_buf+IMG_SZ;

    if (*state)
    {
        img_fir3_t_raw(img_out,img_buf1,img_buf,img_in,coff);
        img_copy(img_buf1,img_in,IMG_SZ);
        *state=0;
    }
    else
    {
        img_fir3_t_raw(img_out,img_buf,img_buf1,img_in,coff);
        img_copy(img_buf,img_in,IMG_SZ);
        *state=1;
    }

    return img_out;
}


/** 
 * @fn              float *img_mid3_t_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2)
 * @details         图像序列的中值滤波,使用前2帧和当前帧数据
 * @param [in]      img0，img1，img2为历史图像帧(指针)，img0对应最老图像，img2对应最新（当前）图像
 * @param [out]     img_out：指针，指向空间存放滤波结果
 * @retval          float *：和img_out相同
 */
float *img_mid3_t_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2)
{
    float *p0=img_in0, *p1=img_in1, *p2=img_in2;
    float *q=img_out,*q_end=img_out+IMG_WID*IMG_HGT;

    for (;q<q_end;p0++,p1++,p2++,q++)
    {
        *q=MID3(*p0,*p1,*p2);
    }
    return img_out;
}


/** 
 * @fn              float *img_mid3_t(float *img_out, float *img_buf, float *img_in, int *state)
 * @details         图像序列的中值滤波,使用前2帧和当前帧数据
 * @param [in]      float *img_buf：历史图像帧(指针)，连续存放最近2帧图像
 * @param [in]      float *img_in：指针，指向最新输入图像
 * @param [out]     float *img_out：指针，指向空间存放滤波结果
 * @param [inout]   int state：指针，指向滤波状态变量，初始值需设为0，指向的内容在运行后被修改
 * @retval          float*：和img_out相同
 */
float *img_mid3_t(float *img_out, float *img_buf, float *img_in, int *state)
{
    float *img_buf1=img_buf+IMG_SZ;

    if (*state)
    {
        img_mid3_t_raw(img_out,img_buf1,img_buf,img_in);
        img_copy(img_buf1,img_in,IMG_SZ);
        *state=0;
    }
    else
    {
        img_mid3_t_raw(img_out,img_buf,img_buf1,img_in);
        img_copy(img_buf,img_in,IMG_SZ);
        *state=1;
    }

    return img_out;
}


/** 
 * @fn              float *img_mid5_t_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2)
 * @details         图像序列的中值滤波,使用前4帧和当前帧数据
 * @param [in]      img0，img1，img2, img3, img4为历史图像帧(指针)，img0对应最老图像，img4对应最新（当前）图像
 * @param [out]     img_out：指针，指向空间存放滤波结果
 * @retval          float *：和img_out相同
 */
float *img_mid5_t_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2, float *img_in3, float *img_in4)
{
    float *p0=img_in0, *p1=img_in1, *p2=img_in2, *p3=img_in4, *p4=img_in4;
    float *q=img_out,*q_end=img_out+IMG_WID*IMG_HGT;

    for (;q<q_end;p0++,p1++,p2++,p3++,p4++,q++)
        *q=mid5(*p0,*p1,*p2,*p3,*p4);
    return img_out;
}


/** 
 * @fn              float *img_mid5_t(float *img_out, float *img_buf, float *img_in, int *state)
 * @details         图像序列的中值滤波,使用前4帧和当前帧数据
 * @param [in]      float *img_buf：历史图像帧(指针)，连续存放最近4帧图像
 * @param [in]      float *img_in：指针，指向最新输入图像
 * @param [out]     float *img_out：指针，指向空间存放滤波结果
 * @param [inout]   int state：指针，指向滤波状态变量(最老的图像帧在img_buf中的位置），初始值需设为0，指向的内容在运行后被修改
 * @retval          float*：和img_out相同
 */
float *img_mid5_t(float *img_out, float *img_buf, float *img_in, int *state)
{
    float *img_buf0=img_buf+IMG_SZ*  (*state);
    float *img_buf1=img_buf+IMG_SZ*(((*state)+1)%4);
    float *img_buf2=img_buf+IMG_SZ*(((*state)+2)%4);
    float *img_buf3=img_buf+IMG_SZ*(((*state)+3)%4);

    img_mid5_t_raw(img_out,img_buf0,img_buf1,img_buf2,img_buf3,img_in);
    img_copy(img_buf0,img_in,IMG_SZ);
    *state=((*state)+1)%4;

    return img_out;
}


/** 
 * @fn              float *img_mid5_avg_t_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2)
 * @details         图像序列的平均中值滤波,使用前4帧和当前帧数据，5帧数据中对应位置像素值，去除最大最小值后平均
 * @param [in]      img0，img1，img2, img3, img4为历史图像帧(指针)，img0对应最老图像，img4对应最新（当前）图像
 * @param [out]     img_out：指针，指向空间存放滤波结果
 * @retval          float *：和img_out相同
 */
float *img_minmax_avg5_t_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2, float *img_in3, float *img_in4)
{
    float *p0=img_in0, *p1=img_in1, *p2=img_in2, *p3=img_in4, *p4=img_in4;
    float *q=img_out,*q_end=img_out+IMG_WID*IMG_HGT;

    for (;q<q_end;p0++,p1++,p2++,p3++,p4++,q++)
        *q=minmax_avg5(*p0,*p1,*p2,*p3,*p4);
    return img_out;
}

/** 
 * @fn              float *img_minmax_avg5_t(float *img_out, float *img_buf, float *img_in, int *state)
 * @details         图像序列的平均中值滤波,使用前4帧和当前帧数据，，5帧数据中对应位置像素，去除最大最小值后平均
 * @param [in]      float *img_buf：历史图像帧(指针)，连续存放最近4帧图像
 * @param [in]      float *img_in：指针，指向最新输入图像
 * @param [out]     float *img_out：指针，指向空间存放滤波结果
 * @param [inout]   int state：指针，指向滤波状态变量(最老的图像帧在img_buf中的位置），初始值需设为0，指向的内容在运行后被修改
 * @retval          float*：和img_out相同
 */
float *img_minmax_avg5_t(float *img_out, float *img_buf, float *img_in, int *state)
{
    float *img_buf0=img_buf+IMG_SZ*  (*state);
    float *img_buf1=img_buf+IMG_SZ*(((*state)+1)%4);
    float *img_buf2=img_buf+IMG_SZ*(((*state)+2)%4);
    float *img_buf3=img_buf+IMG_SZ*(((*state)+3)%4);

    img_minmax_avg5_t_raw(img_out,img_buf0,img_buf1,img_buf2,img_buf3,img_in);
    img_copy(img_buf0,img_in,IMG_SZ);
    *state=((*state)+1)%4;

    return img_out;
}


/** 
 * @fn              float *img_mid_cross(float *img_out, float *img_in)
 * @details         图像空间域中值滤波，使用十字滤波模板（共5个像素）。注意：滤波输出图像的最外圈边沿（1层像素）是无效数据。滤波器模板如下：
 *                             0(p)
 *                   1(p+W-1)  2(p+W)  3(p+W+1)
 *                             4(p+2W)
 *                  注意：滤波输出图像的最外圈边沿（1层像素）是无效数据
 * @param [in]      float *img_in：指针，指向待滤波图像
 * @param [out]     float *img_out：指针，指向的空间存放图像运算结果
 * @retval          float *：和img_out相同
 */ 
float *img_mid_cross(float *img_out, float *img_in)
{
    float *p0=img_in+1;
    float *p1=p0+IMG_WID-1;
    float *p2=p0+IMG_WID  ;
    float *p3=p0+IMG_WID+1;
    float *p4=p0+IMG_WID+IMG_WID;
    
    float *q=img_out+IMG_WID+1,*q_end=img_out+IMG_WID*(IMG_HGT-1)-1;
    
    for (;q<q_end;q++,p0++,p1++,p2++,p3++,p4++)
        *q=mid5(*p0,*p1,*p2,*p3,*p4);

    return img_out;
}


/** 
 * @fn              float *img_mid_cross_sa(float *img_inout, float *img_in)
 * @brief           图像空间域中值滤波（原址运算），使用十字滤波模板（共5个像素）。注意：滤波输出图像的最外圈边沿（1层像素）是无效数据
 * @details         滤波器模板如下
 *                             0(p)
 *                   1(p+W-1)  2(p+W)  3(p+W+1)
 *                             4(p+2W)
 *                  注意：滤波输出图像的最外圈边沿（1层像素）是无效数据
 * @param [in]      float *img_inout：指针，指向待滤波图像和图像运算结果
 * @param [inout]   ring_buf_f32_s *rbuf：环形缓冲器，存放3行数据
 * @param [inout]   float *img_inout：指针，指向的空间存放图像运算结果
 * @retval          float *：和img_inout相同
 */ 
float *img_mid_cross_sa(float *img_inout, float *img_in, struct ring_buf_f32_s *rbuf)
{
    float *p0=img_in+1;
    float *p1=p0+IMG_WID-1;
    float *p2=p0+IMG_WID  ;
    float *p3=p0+IMG_WID+1;
    float *p4=p0+IMG_WID+IMG_WID;
    
    float *q =img_inout+IMG_WID+1-3*IMG_WID,*q_end=img_inout+IMG_WID*(IMG_HGT-1)-1-3*IMG_WID;

    int n=3*IMG_WID;
    
    float s;
    for (;q<q_end;q++,p0++,p1++,p2++,p3++,p4++)
    {
        s=mid5(*p0,*p1,*p2,*p3,*p4);
        if (n)
        {
            n--;
            ring_buf_f32_io(rbuf,s);
        }
        else
            *q=ring_buf_f32_io(rbuf,s);
    }
    for (n=0;n<3*IMG_WID;n++,q++)
        *q=ring_buf_f32_io(rbuf,0);
    return img_inout;
}


/** 
 * @fn              float *img_iir_sos(float *img_out, float *img_in, float img_st0, float *img_st1, float *coff)
 * @details         使用2阶IIR滤波器的图像时域滤波，
 *                  Maltab的SOS矩阵数据格式是（每行数据格式） [b1 b2 b3, a1 a2 a3], 由于a1=1，因此在填入下面的系数数组时被去除
 *                  Maltab的G里面是sc（尺度缩放）数据
 * @param [in]      float *img_in：指针，指向待滤波图像
 * @param [inout]   float *img_st1,*img_st2：指针，指向滤波状态数据（图像）
 * @param [in]      float *coff：指针，指向滤波加权系数数组{b1,b2,b3,a2,a3,sc}
 * @param [out]     float *img_out：指针，指向的空间存放图像运算结果
 * @retval          float *：和img_inout相同
 */ 
float *img_iir_sos(float *img_out, float *img_in, float *img_st1, float *img_st2, float *coff)
{
    float b1=*coff,b2=*(coff+1),b3=*(coff+2),a2=*(coff+3),a3=*(coff+4),sc=*(coff+5);

    // img_out[:]=b1*img_in+img_st1[:]
    img_mul_f32(img_out,img_in,b1);
    img_cum(img_out,img_st1);

    // img_st1[:]=b1*img_in[:]+img_st2[:]-a2*img_out[:]
    img_mul_f32(img_st1,img_in,b2);
    img_cum(img_st1,img_st2);
    img_mac(img_st1,img_out,-a2);

    // img_st2[:]=b3*img_in[:]-a3*img_out[:]
    img_mul_f32(img_st2,img_in,b3);
    img_mac(img_st2,img_out,-a3);

    // int_out[:]*=sc
    img_prod_f32(img_out,sc);

    return img_out;
};


/** 
 * @fn              float *img_weighted_iir(float *img_out, float *img_in, float *img_in_w_avg, float *img_w, float *img_w_avg, float alpha)
 * @details         图像加权IIR平均，使用以下算法:
 *                  img_w_avg[:]=img_w_avg[:]*alpha+(1-alpha)img_w[:]
 *                  img_in_w_avg[:]=img_in_w_avg[:]*alpha+(1-alpha)img_w[:].*img_in[:]
 *                  img_out[:]=img_in_w_avg[:]./img_w_avg[:]
 * @param [in]      float *img_in：指针，指向待滤波图像
 * @param [in]      float *img_w：指针，指向加权数据
 * @param [inout]   float *img_in_w_avg：指针，指向滤波状态数据,内容在该函数运行后更新
 * @param [inout]   float *img_w_avg：指针，指向滤波状态数据,内容在该函数运行后更新
 * @param [in]      float alpha：滤波器遗忘因子(0~1)越接近0，越“健忘”
 * @param [out]     float *img_out：指针，指向的空间存放图像运算结果
 * @retval          float *：和img_inout相同
 */ 
float *img_weighted_iir(float *img_out, float *img_in, float *img_in_w_avg, float *img_w, float *img_w_avg, float alpha)
{
    float *q=img_out, *q_end=img_out+IMG_SZ;

    float *p1=img_in;
    float *p2=img_in_w_avg;
    float *p3=img_w;
    float *p4=img_w_avg;

    for (;q<q_end;q++,p1++,p2++,p3++,p4++)
    {
        *p4=(*p4)*alpha+(float)(1.0-alpha)*(*p3);
        *p2=(*p2)*alpha+(float)(1.0-alpha)*(*p3)*(*p1);
        if (*p4)
            *q=(*p2)/(*p4);
        else
            *q=0;
    }

    return img_out;
}


/** 
 * @fn              float *img_max3_t_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2)
 * @details         从3帧图像序列中找到的每个位置的像素最大值,使用前2帧和当前帧数据
 * @param [in]      float* img0，img1，img2为历史图像帧(指针)，img0对应最老图像，img2对应最新（当前）图像
 * @param [out]     float* img_out：指针，指向空间存放滤波结果
 * @retval          float *：和img_out相同
 */
float *img_max3_t_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2)
{
    float *p0=img_in0, *p1=img_in1, *p2=img_in2;
    float *q=img_out,*q_end=img_out+IMG_WID*IMG_HGT;

    for (;q<q_end;p0++,p1++,p2++,q++)
    {
        *q=MAX3(*p0,*p1,*p2);
    }
    return img_out;
}


/** 
 * @fn              float *img_max3_t(float *img_out, float *img_buf, float *img_in, int *state)
 * @details         从连续输入的最近3帧图像序列中找到的每个位置的像素最大值,使用前2帧和当前帧数据
 * @param [in]      float *img_buf：历史图像帧(指针)，连续存放最近2帧图像
 * @param [in]      float *img_in：指针，指向最新输入图像
 * @param [out]     float *img_out：指针，指向空间存放滤波结果
 * @param [inout]   int state：指针，指向滤波状态变量，初始值需设为0，指向的内容在运行后被修改
 * @retval          float*：和img_out相同
 */
float *img_max3_t(float *img_out, float *img_buf, float *img_in, int *state)
{
    float *img_buf1=img_buf+IMG_SZ;

    if (*state)
    {
        img_max3_t_raw(img_out,img_buf1,img_buf,img_in);
        img_copy(img_buf1,img_in,IMG_SZ);
        *state=0;
    }
    else
    {
        img_max3_t_raw(img_out,img_buf,img_buf1,img_in);
        img_copy(img_buf,img_in,IMG_SZ);
        *state=1;
    }

    return img_out;
}


/** 
 * @fn              float *img_max5_t_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2)
 * @details         从连续输入的最近5帧图像序列中找到的每个位置的像素最大值,使用前4帧和当前帧数据
 * @param [in]      float *img0，img1，img2, img3, img4为历史图像帧(指针)，img0对应最老图像，img4对应最新（当前）图像
 * @param [out]     float *img_out：指针，指向空间存放滤波结果
 * @retval          float *：和img_out相同
 */
float *img_max5_t_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2, float *img_in3, float *img_in4)
{
    float *p0=img_in0, *p1=img_in1, *p2=img_in2, *p3=img_in4, *p4=img_in4;
    float *q=img_out,*q_end=img_out+IMG_WID*IMG_HGT;

    for (;q<q_end;p0++,p1++,p2++,p3++,p4++,q++)
        *q=max5(*p0,*p1,*p2,*p3,*p4);
    return img_out;
}


/** 
 * @fn              float *img_max5_t(float *img_out, float *img_buf, float *img_in, int *state)
 * @details         从连续输入的最近5帧图像序列中找到的每个位置的像素最大值,使用前4帧和当前帧数据
 * @param [in]      float *img_buf：历史图像帧(指针)，连续存放最近4帧图像
 * @param [in]      float *img_in：指针，指向最新输入图像
 * @param [out]     float *img_out：指针，指向空间存放滤波结果
 * @param [inout]   int state：指针，指向滤波状态变量(最老的图像帧在img_buf中的位置），初始值需设为0，指向的内容在运行后被修改
 * @retval          float*：和img_out相同
 */
float *img_max5_t(float *img_out, float *img_buf, float *img_in, int *state)
{
    float *img_buf0=img_buf+IMG_SZ*  (*state);
    float *img_buf1=img_buf+IMG_SZ*(((*state)+1)%4);
    float *img_buf2=img_buf+IMG_SZ*(((*state)+2)%4);
    float *img_buf3=img_buf+IMG_SZ*(((*state)+3)%4);

    img_max5_t_raw(img_out,img_buf0,img_buf1,img_buf2,img_buf3,img_in);
    img_copy(img_buf0,img_in,IMG_SZ);
    *state=((*state)+1)%4;

    return img_out;
}


/** 
 * @fn              float *img_min5_t_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2)
 * @details         从连续输入的最近5帧图像序列中找到的每个位置的像素最小值,使用前4帧和当前帧数据
 * @param [in]      float *img0，img1，img2, img3, img4为历史图像帧(指针)，img0对应最老图像，img4对应最新（当前）图像
 * @param [out]     float *img_out：指针，指向空间存放滤波结果
 * @retval          float *：和img_out相同
 */
float *img_min5_t_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2, float *img_in3, float *img_in4)
{
    float *p0=img_in0, *p1=img_in1, *p2=img_in2, *p3=img_in4, *p4=img_in4;
    float *q=img_out,*q_end=img_out+IMG_WID*IMG_HGT;

    for (;q<q_end;p0++,p1++,p2++,p3++,p4++,q++)
        *q=max5(*p0,*p1,*p2,*p3,*p4);
    return img_out;
}


/** 
 * @fn              float *img_min5_t(float *img_out, float *img_buf, float *img_in, int *state)
 * @details         从连续输入的最近4帧图像序列中找到的每个位置的像素最小值,使用前4帧和当前帧数据
 * @param [in]      float *img_buf：历史图像帧(指针)，连续存放最近4帧图像
 * @param [in]      float *img_in：指针，指向最新输入图像
 * @param [out]     float *img_out：指针，指向空间存放滤波结果
 * @param [inout]   int state：指针，指向滤波状态变量(最老的图像帧在img_buf中的位置），初始值需设为0，指向的内容在运行后被修改
 * @retval          float*：和img_out相同
 */
float *img_min5_t(float *img_out, float *img_buf, float *img_in, int *state)
{
    float *img_buf0=img_buf+IMG_SZ*  (*state);
    float *img_buf1=img_buf+IMG_SZ*(((*state)+1)%4);
    float *img_buf2=img_buf+IMG_SZ*(((*state)+2)%4);
    float *img_buf3=img_buf+IMG_SZ*(((*state)+3)%4);

    img_max5_t_raw(img_out,img_buf0,img_buf1,img_buf2,img_buf3,img_in);
    img_copy(img_buf0,img_in,IMG_SZ);
    *state=((*state)+1)%4;

    return img_out;
}


/** 
 * @fn              float *img_nnf_sqr3(float *img_out, float *img_in, float *coff)
 * @detai           像素最近邻选择滤波，使用使用3x3滤波模板，如果邻近像素和中心像素差异超过门限则使用空间十字模板（5点）中值滤波
 *                  最近邻像素的位置为3x3矩阵，如下所示：
 *                  0(p)    1(p+1)    2(p+2)
 *                  3(p+W)  4(p+W+1)  5(p+W+2)
 *                  6(P+2W) 7(p+2W+1) 8(p+2W+2)
 *                  注意：滤波输出图像的最外圈边沿（1层像素）是无效数据
 *                  计算步骤为：1). 计算3x3邻近像素差别；2). 对于超过门限的点，用十字模板（5个点）的中值取代
 * @param [in]      float *img_in：指针，指向待滤波图像
 * @param [in]      float th：滤波门限
 * @param [out]     float *img_out：指针，指向的空间存放图像运算结果
 * @retval          float *：和img_out相同
 */ 
float *img_nnf_sqr3(float *img_out, float *img_in, float th)
{
    float s;
    float *p0=img_in          , *p1=img_in          +1, *p2=img_in          +2;
    float *p3=img_in+  IMG_WID, *p4=img_in+  IMG_WID+1, *p5=img_in+  IMG_WID+2;
    float *p6=img_in+2*IMG_WID, *p7=img_in+2*IMG_WID+1, *p8=img_in+2*IMG_WID+2;

    float *q=img_out+IMG_WID+1,*q_end=img_out+IMG_WID*(IMG_HGT-1)-1;

    for (;q<q_end;q++,p0++,p1++,p2++,p3++,p4++,p5++,p6++,p7++,p8++)
    {
        s=min8((float)fabs((*p0)-(*p4)),(float)fabs((*p1)-(*p4)),(float)fabs((*p2)-(*p4)),(float)fabs((*p3)-(*p4)),
               (float)fabs((*p5)-(*p4)),(float)fabs((*p6)-(*p4)),(float)fabs((*p7)-(*p4)),(float)fabs((*p8)-(*p4)));
        if (s<th)
            *q=*p4;
        else
            *q=mid5(*p1,*p3,*p4,*p5,*p7);
    }

    return img_out;
}


/** 
 * @fn              float *img_fb_mid5_t_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2)
 * @details         图像序列的前向后向选择中值滤波,使用前4帧和当前帧数据
 *                  步骤为：1)计算包括当前帧的前3帧中值（前中值），和包括当前帧的后3帧中值（后中值）; 
 *                  2)计算两个中值的差，超过门限时，用后中值取代当前点，否则用5帧（前后各2帧加上当前帧）的5个点中值代替当前帧像素         
 * @param [in]      float *img0，img1，img2, img3, img4为历史图像帧(指针)，img0对应最老图像，img4对应最新（当前）图像
 * @param [in]      float th使用前向MID3滤波结果（新数据）的门限
 * @param [out]     float *img_out：指针，指向空间存放滤波结果
 * @retval          float *：和img_out相同
 */
float *img_fb_mid3_t_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2, float *img_in3, float *img_in4, float th)
{
    float *p0=img_in0, *p1=img_in1, *p2=img_in2, *p3=img_in4, *p4=img_in4;
    float *q=img_out,*q_end=img_out+IMG_WID*IMG_HGT;

    float a,b;
    for (;q<q_end;p0++,p1++,p2++,p3++,p4++,q++)
    {
        a=MID3(*p0,*p1,*p2);
        b=MID3(*p2,*p3,*p4);
        if (fabs(a-b)>th)
            *q=*p2;
        else
            *q=mid5(*p0,*p1,*p2,*p3,*p4);
    }
    return img_out;
}


/** 
 * @fn              float *img_fb_mid5_t(float *img_out, float *img_buf, float *img_in, int *state)
 * @details         图像序列的前向后向选择中值滤波,使用前4帧和当前帧数据
 *                  步骤为：1)计算包括当前帧的前3帧中值（前中值），和包括当前帧的后3帧中值（后中值）; 
 *                  2)计算两个中值的差，超过门限时，用后中值取代当前点，否则用5帧（前后各2帧加上当前帧）的5个点中值代替当前帧像素         
 * @param [in]      float *img0，img1，img2, img3, img4为历史图像帧(指针)，img0对应最老图像，img4对应最新（当前）图像
 * @param [in]      float *img_buf：历史图像帧(指针)，连续存放最近4帧图像
 * @param [in]      float *img_in：指针，指向最新输入图像
 * @param [in]      float th使用前向MID3滤波结果（新数据）的门限
 * @param [out]     float *img_out：指针，指向空间存放滤波结果
 * @param [inout]   int state：指针，指向滤波状态变量(最老的图像帧在img_buf中的位置），初始值需设为0，指向的内容在运行后被修改
 * @retval          float*：和img_out相同
 */
float *img_fb_mid3_t(float *img_out, float *img_buf, float *img_in, float th, int *state)
{
    float *img_buf0=img_buf+IMG_SZ*  (*state);
    float *img_buf1=img_buf+IMG_SZ*(((*state)+1)%4);
    float *img_buf2=img_buf+IMG_SZ*(((*state)+2)%4);
    float *img_buf3=img_buf+IMG_SZ*(((*state)+3)%4);

    img_fb_mid3_t_raw(img_out,img_buf0,img_buf1,img_buf2,img_buf3,img_in,th);
    img_copy(img_buf0,img_in,IMG_SZ);
    *state=((*state)+1)%4;

    return img_out;
}


/** 
 * @fn              float *img_nnf_sqr3_mid5_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2, float *img_in3, float *img_in4, float th)
 * @details         像素最近邻选择滤波，使用使用3x3滤波模板，如果邻近像素和中心像素差异过大则使用5帧图像的时间中值滤波
 *                  滤波器模板如下
 *                  0(p)    1(p+1)    2(p+2)
 *                  3(p+W)  4(p+W+1)  5(p+W+2)
 *                  6(P+2W) 7(p+2W+1) 8(p+2W+2)
 *                  注意：滤波输出图像的最外圈边沿（1层像素）是无效数据
 *                  步骤为：1). 计算计算当前像素和周围3x3邻近像素差别;
 *                  2). 对于超过门限的点，使用5帧图像的时间中值滤波
 * @param [in]      float *img_in：指针，指向待滤波图像
 * @param [in]      float th：滤波门限
 * @param [out]     float *img_out：指针，指向的空间存放图像运算结果
 * @retval          float *：和img_out相同
 */ 
float *img_nnf_sqr3_mid5_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2, float *img_in3, float *img_in4, float th)
{
    float s;
    float *p0=img_in2          , *p1=img_in2          +1, *p2=img_in2          +2;
    float *p3=img_in2+  IMG_WID, *p4=img_in2+  IMG_WID+1, *p5=img_in2+  IMG_WID+2;
    float *p6=img_in2+2*IMG_WID, *p7=img_in2+2*IMG_WID+1, *p8=img_in2+2*IMG_WID+2;

    float *r0=img_in0+IMG_WID+1, *r1=img_in1+IMG_WID+1, *r2=img_in2+IMG_WID+1, *r3=img_in3+IMG_WID+1, *r4=img_in4+IMG_WID+1;

    float *q=img_out+IMG_WID+1,*q_end=img_out+IMG_WID*(IMG_HGT-1)-1;

    for (;q<q_end;q++,p0++,p1++,p2++,p3++,p4++,p5++,p6++,p7++,p8++,r0++,r1++,r2++,r3++,r4++)
    {
        s=min8((float)fabs((*p0)-(*p4)),(float)fabs((*p1)-(*p4)),(float)fabs((*p2)-(*p4)),(float)fabs((*p3)-(*p4)),
               (float)fabs((*p5)-(*p4)),(float)fabs((*p6)-(*p4)),(float)fabs((*p7)-(*p4)),(float)fabs((*p8)-(*p4)));
        if (s<th)
            *q=*p4;
        else
            *q=mid5(*r0,*r1,*r2,*r3,*r4);
    }

    return img_out;
}


/** 
 * @fn              float *img_nnf_sqr3_mid5(float *img_out, float *img_buf, float *img_in, int *state,float th)
 * @details         像素最近邻选择滤波，使用使用3x3滤波模板，如果邻近像素和中心像素差异过大则使用5帧图像的时间中值滤波
 *                  步骤为：1). 计算当前像素和周围3x3邻近像素差别; 
 *                  2). 对于超过门限的点，使用5帧图像的时间中值滤波
 * @param [in]      float *img_buf：历史图像帧(指针)，连续存放最近4帧图像
 * @param [in]      float *img_in：指针，指向最新输入图像
 * @param [in]      float th：滤波门限
 * @param [out]     float *img_out：指针，指向空间存放滤波结果
 * @param [inout]   int state：指针，指向滤波状态变量(最老的图像帧在img_buf中的位置），初始值需设为0，指向的内容在运行后被修改
 * @retval          float*：和img_out相同
 */
float *img_nnf_sqr3_mid5(float *img_out, float *img_buf, float *img_in, int *state,float th)
{
    float *img_buf0=img_buf+IMG_SZ*  (*state);
    float *img_buf1=img_buf+IMG_SZ*(((*state)+1)%4);
    float *img_buf2=img_buf+IMG_SZ*(((*state)+2)%4);
    float *img_buf3=img_buf+IMG_SZ*(((*state)+3)%4);

    img_nnf_sqr3_mid5_raw(img_out,img_buf0,img_buf1,img_buf2,img_buf3,img_in,th);
    img_copy(img_buf0,img_in,IMG_SZ);
    *state=((*state)+1)%4;

    return img_out;
}


/** 
 * @fn              float *img_mid7_st_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2)
 * @details         图像序列的时空中值滤波,使用2帧历史数据
 *                  从前后帧和当前帧得到7个像素，用中值取代当前帧数据。当前帧像素点位置如下
 *                            0(p)
 *                  1(p+W-1)  2(p+W)  3(p+W+1)
 *                            4(p+2W)
 *                  前后一帧使用2号位置像素数数据，共7个像素数据
 * @param [in]      img0，img1，img2为历史图像帧(指针)，img0对应最老图像，img2对应最新（当前）图像
 * @param [out]     img_out：指针，指向空间存放滤波结果
 * @retval          float *：和img_out相同
 */
float *img_mid7_st_raw(float *img_out, float *img_in0, float *img_in1, float *img_in2)
{
    // 当前图
    float *p0=img_in1+1;
    float *p1=img_in1+IMG_WID;
    float *p2=img_in1+IMG_WID+1; // 中心点
    float *p3=img_in1+IMG_WID+2;
    float *p4=img_in1+2*IMG_WID+1;
    
    // 前后图
    float *p5=img_in0+IMG_WID+1;    // 前图中心点
    float *p6=img_in2+IMG_WID+1;    // 后图中心点

    // 输出指针
    float *q=img_out+IMG_WID+1,*q_end=img_out+IMG_WID*(IMG_HGT-1)-1;
    
    //中值滤波
    for (;q<q_end;q++,p0++,p1++,p2++,p3++,p4++,p5++,p6++)
        *q=mid7(*p0,*p1,*p2,*p3,*p4,*p5,*p5);

    return img_out;
}


/** 
 * @fn              float *img_mid7_st(float *img_out, float *img_buf, float *img_in, int *state)
 * @details         从前后帧和当前帧得到7个像素，用中值取代当前帧数据。当前帧像素点位置如下
 *                            0(p)
 *                  1(p+W-1)  2(p+W)  3(p+W+1)
 *                            4(p+2W)
 *                  前后一帧使用2号位置像素数数据，共7个像素数据
 * @param [in]      float *img_buf：历史图像帧(指针)，连续存放最近2帧图像
 * @param [in]      float *img_in：指针，指向最新输入图像
 * @param [out]     float *img_out：指针，指向空间存放滤波结果
 * @param [inout]   int state：指针，指向滤波状态变量，初始值需设为0，指向的内容在运行后被修改
 * @retval          float*：和img_out相同
 */
float *img_mid7_st(float *img_out, float *img_buf, float *img_in, int *state)
{
    float *img_buf1=img_buf+IMG_SZ;

    if (*state)
    {
        img_mid7_st_raw(img_out,img_buf1,img_buf,img_in);
        img_copy(img_buf1,img_in,IMG_SZ);
        *state=0;
    }
    else
    {
        img_mid7_st_raw(img_out,img_buf,img_buf1,img_in);
        img_copy(img_buf,img_in,IMG_SZ);
        *state=1;
    }

    return img_out;
}


static float sqr_f32(float x) { return x*x; }

// 平面匹配滤波器，9点, 
// 输入：3x3=9个像素点深度，
//      z0 z1 z2
//      z3 z4 z5
//      z6 z7 z8
// 计算原理：
//    从9个点中，找出6个点，计算拟合的平面离那6和点的距离误差，找出最匹配的6个点，作为匹配结果，修正中间点(z4)的深度 
float img_plane_mf_pix(float z0, float z1, float z2, float z3, float z4, float z5, float z6, float z7, float z8)
{
    float e0,e1,e2,e3,e4,e5,e6,e7;
    float minv;
    int min_id;
    float zc;

    e0= (sqr_f32(-5*z0+4*z1+  z2+3*z3     -3*z5)+   //  z0 z1 z2    // 用上方6个点你和平面，计算拟合误差e0
         sqr_f32( 4*z0-8*z1+4*z2               )+   //  z3 z4 z5
         sqr_f32(   z0+4*z1-5*z2-3*z3     +3*z5)+   //  *  *  *
         sqr_f32( 3*z0     -3*z2-5*z3+4*z4+  z5)+       
         sqr_f32(                4*z3-8*z4+4*z5)+       
         sqr_f32(-3*z0     +3*z2+  z3+4*z4-5*z5))/144;

    e1= (sqr_f32(-3*z0+3*z1-  z2+3*z4-  z5-  z8)+   //  z0 z1 z2    // 用右上方6个点你和平面，计算拟合误差e1
         sqr_f32( 3*z0-7*z1+3*z2+  z4+  z5-  z8)+   //  *  z4 z5
         sqr_f32(  -z0+3*z1-3*z2-  z4+3*z5-  z8)+   //  *  *  z8
         sqr_f32( 3*z0+  z1-  z2-7*z4+  z5+3*z8)+
         sqr_f32(  -z0+  z1+3*z2+  z4-7*z5+3*z8)+
         sqr_f32(  -z0-  z1-  z2+3*z4+3*z5-3*z8))/100;
    
    e2= (sqr_f32(-5*z1+3*z2+4*z4     +  z7-3*z8)+   //  *  z1 z2    // 用右方6个点你和平面，计算拟合误差e2
         sqr_f32( 3*z1-5*z2     +4*z5-3*z7+  z8)+   //  *  z4 z5
         sqr_f32( 4*z1     -8*z4     +4*z7     )+   //  *  z7 z8
         sqr_f32(      4*z2     -8*z5     +4*z8)+
         sqr_f32(   z1-3*z2+4*z4     -5*z7+3*z8)+
         sqr_f32(-3*z1+  z2     +4*z5+3*z7-5*z8))/144;

    e3= (sqr_f32(-3*z2+3*z4+3*z5-  z6-  z7-  z8)+   //  *  *  z2    // 用右下方6个点你和平面，计算拟合误差e3
         sqr_f32( 3*z2-7*z4+  z5+3*z6+  z7-  z8)+   //  *  z4 z5
         sqr_f32( 3*z2+  z4-7*z5-  z6+  z7+3*z8)+   //  z6 z7 z8
         sqr_f32(  -z2+3*z4-  z5-3*z6+3*z7-  z8)+
         sqr_f32(  -z2+  z4+  z5+3*z6-7*z7+3*z8)+
         sqr_f32(  -z2-  z4+3*z5-  z6+3*z7-3*z8))/100;

    e4= (sqr_f32(-5*z3+4*z4+  z5+3*z6     -3*z8)+   //  *  *  *     // 用下方6个点你和平面，计算拟合误差e4
         sqr_f32( 4*z3-8*z4+4*z5               )+   //  z3 z4 z5
         sqr_f32(   z3+4*z4-5*z5-3*z6     +3*z8)+   //  z6 z7 z8
         sqr_f32( 3*z3     -3*z5-5*z6+4*z7+  z8)+
         sqr_f32(                4*z6-8*z7+4*z8)+
         sqr_f32(-3*z3     +3*z5+  z6+4*z7-5*z8))/144;

    e5= (sqr_f32(-3*z0+3*z3+3*z4-  z6-  z7-  z8)+   //  z0 *  *     // 用左下方6个点你和平面，计算拟合误差e5
         sqr_f32( 3*z0-7*z3+  z4+3*z6+  z7-  z8)+   //  z3 z4 *
         sqr_f32( 3*z0+  z3-7*z4-  z6+  z7+3*z8)+   //  z6 z7 z8
         sqr_f32(  -z0+3*z3-  z4-3*z6+3*z7-  z8)+
         sqr_f32(  -z0+  z3+  z4+3*z6-7*z7+3*z8)+
         sqr_f32(  -z0-  z3+3*z4-  z6+3*z7-3*z8))/100;

    e6= (sqr_f32(-5*z0+3*z1+4*z3     +  z6-3*z7)+   //  z0 z2 *     // 用左方6个点你和平面，计算拟合误差e6
         sqr_f32( 3*z0-5*z1     +4*z4-3*z6+  z7)+   //  z3 z4 *
         sqr_f32( 4*z0     -8*z3     +4*z6     )+   //  z6 z7 *
         sqr_f32(      4*z1     -8*z4     +4*z7)+
         sqr_f32(   z0-3*z1+4*z3     -5*z6+3*z7)+
         sqr_f32(-3*z0+  z1     +4*z4+3*z6-5*z7))/144;

    e7= (sqr_f32(-3*z0+3*z1-  z2+3*z3-  z4-  z6)+   //  z0 z1 z2    // 用左上6个点你和平面，计算拟合误差e7
         sqr_f32( 3*z0-7*z1+3*z2+  z3+  z4-  z6)+   //  z3 z4 *
         sqr_f32(-  z0+3*z1-3*z2-  z3+3*z4-  z6)+   //  z6 *  * 
         sqr_f32( 3*z0+  z1-  z2-7*z3+  z4+3*z6)+
         sqr_f32(  -z0+  z1+3*z2+  z3-7*z4+3*z6)+
         sqr_f32(  -z0-  z1-  z2+3*z3+3*z4-3*z6))/100;

    // 根据e0~e7，找到最优拟合方案
    minv=e0;
    min_id=0;
    if (e1<minv) { minv=e1; min_id=1; }
    if (e2<minv) { minv=e2; min_id=2; }
    if (e3<minv) { minv=e3; min_id=3; }
    if (e4<minv) { minv=e4; min_id=4; }
    if (e5<minv) { minv=e5; min_id=5; }
    if (e6<minv) { minv=e6; min_id=6; }
    if (e7<minv) { minv=e7; min_id=7; }

    // 计算中心点深度修正结果
         if (min_id==0) zc=(                  z3+  z4+  z5)/3 ;
    else if (min_id==1) zc=( 3*z0+  z1-  z2+3*z4+  z5+3*z8)/10;
    else if (min_id==2) zc=(   z1     +  z4     +  z7     )/3 ;    
    else if (min_id==3) zc=( 3*z2+3*z4+  z5+3*z6+  z7  -z8)/10;    
    else if (min_id==4) zc=(   z3+  z4+  z5               )/3 ;    
    else if (min_id==5) zc=( 3*z0+  z3+3*z4  -z6+  z7+3*z8)/10;    
    else if (min_id==6) zc=(        z1     +  z4     +  z7)/3 ;    
    else if (min_id==7) zc=(  -z0+  z1+3*z2+  z3+3*z4+3*z6)/10;    

    return zc;
}


float *img_plane_mf_sqr3(float *img_out, float *img_in)
{
    float *p0=img_in          , *p1=img_in          +1, *p2=img_in          +2;
    float *p3=img_in+  IMG_WID, *p4=img_in+  IMG_WID+1, *p5=img_in+  IMG_WID+2;
    float *p6=img_in+2*IMG_WID, *p7=img_in+2*IMG_WID+1, *p8=img_in+2*IMG_WID+2;

    float *q=img_out+IMG_WID+1,*q_end=img_out+IMG_WID*(IMG_HGT-1)-1;

    for (;q<q_end;q++,p0++,p1++,p2++,p3++,p4++,p5++,p6++,p7++,p8++)
        *q=img_plane_mf_pix(*p0,*p1,*p2,*p3,*p4,*p5,*p6,*p7,*p8);

    return img_out;
}

float *img_plane_mf_sqr3_sa(float *img_inout, struct ring_buf_f32_s *rbuf)
{
    float s;
    float *p0=img_inout          , *p1=img_inout          +1, *p2=img_inout          +2;
    float *p3=img_inout+  IMG_WID, *p4=img_inout+  IMG_WID+1, *p5=img_inout+  IMG_WID+2;
    float *p6=img_inout+2*IMG_WID, *p7=img_inout+2*IMG_WID+1, *p8=img_inout+2*IMG_WID+2;

    float *q=img_inout+IMG_WID+1-3*IMG_WID,*q_end=img_inout+IMG_WID*(IMG_HGT-1)-1-3*IMG_WID;
    
    int n=3*IMG_WID;

    for (;q<q_end;q++,p0++,p1++,p2++,p3++,p4++,p5++,p6++,p7++,p8++)
    {
        s=img_plane_mf_pix(*p0,*p1,*p2,*p3,*p4,*p5,*p6,*p7,*p8);

        if (n)
        {
            n--;
            ring_buf_f32_io(rbuf,s);
        }
        else
            *q=ring_buf_f32_io(rbuf,s);
    }

    for (n=0;n<3*IMG_WID;n++,q++)
        *q=ring_buf_f32_io(rbuf,0);

    return img_inout;
}


/** 
 * @fn              float *img_hole_fill(float *img_out, float *img_in, uint8_t *img_mask);
 * @details         像素空洞检测滤波，如果某个像素无效，且他的邻近像素超过(包括）5个非零，则用有效像素平均值填充
 *                  3x3图像滤波器模板如下
 *                  0(p)    1(p+1)    2(p+2)
 *                  3(p+W)  4(p+W+1)  5(p+W+2)
 *                  6(P+2W) 7(p+2W+1) 8(p+2W+2)
 *                  注意：滤波输出图像的最外圈边沿（1层像素）是无效数据
 * @param [in]      float *img_in：指针，指向待滤波图像
 * @param [in]      float *coff：指针，指向9个滤波系数
 * @param [out]     float *img_out：指针，指向的空间存放图像运算结果
 * @param [inout]   uint8_t *img_mask：指针，指向的空间存放空洞指示，注意，填补空洞后会修改该指针对应空间内容
 * @retval          float *：和img_out相同
 */ 
float *img_hole_fill(float *img_out, float *img_in, uint8_t *img_mask)
{
    float *p0=img_in          , *p1=img_in          +1, *p2=img_in          +2;
    float *p3=img_in+  IMG_WID, *p4=img_in+  IMG_WID+1, *p5=img_in+  IMG_WID+2;
    float *p6=img_in+2*IMG_WID, *p7=img_in+2*IMG_WID+1, *p8=img_in+2*IMG_WID+2;

    uint8_t *r0=img_mask          , *r1=img_mask          +1, *r2=img_mask          +2;
    uint8_t *r3=img_mask+  IMG_WID, *r4=img_mask+  IMG_WID+1, *r5=img_mask+  IMG_WID+2;
    uint8_t *r6=img_mask+2*IMG_WID, *r7=img_mask+2*IMG_WID+1, *r8=img_mask+2*IMG_WID+2;

    float *q=img_out+IMG_WID+1,*q_end=img_out+IMG_WID*(IMG_HGT-1)-1;
    
    int k=0;
    
    img_copy(img_out,img_in,IMG_SZ);

    for (;q<q_end;q++,p0++,p1++,p2++,p3++,p4++,p5++,p6++,p7++,p8++,
                      r0++,r1++,r2++,r3++,r4++,r5++,r6++,r7++,r8++)
    {
        if (*r4) continue;  // 非空洞
        
        // 空洞处理
        // 计算空洞的临近像素有效率
        k=*r0+*r1+*r2+*r3+*r5+*r6+*r7+*r8;
        if (k>5)
        {
            *q=0;
            if (*r0) *q+=*p0;
            if (*r1) *q+=*p1;
            if (*r2) *q+=*p2;
            if (*r3) *q+=*p3;
            if (*r5) *q+=*p5;
            if (*r6) *q+=*p6;
            if (*r7) *q+=*p7;
            if (*r8) *q+=*p8;
            *q/=(float)k;
            *r4=1;
        }
    }
    return img_out;
}


// 计算和周围3x3领域点的像素值差的（绝对值）最小值
float *img_nnd_sqr3(float *img_out, float *img_in)
{
    float *p0=img_in          , *p1=img_in          +1, *p2=img_in          +2;
    float *p3=img_in+  IMG_WID, *p4=img_in+  IMG_WID+1, *p5=img_in+  IMG_WID+2;
    float *p6=img_in+2*IMG_WID, *p7=img_in+2*IMG_WID+1, *p8=img_in+2*IMG_WID+2;

    float *q=img_out+IMG_WID+1,*q_end=img_out+IMG_WID*(IMG_HGT-1)-1;

    for (;q<q_end;q++,p0++,p1++,p2++,p3++,p4++,p5++,p6++,p7++,p8++)
        *q=min8((float)fabs((*p0)-(*p4)),(float)fabs((*p1)-(*p4)),(float)fabs((*p2)-(*p4)),(float)fabs((*p3)-(*p4)),
                (float)fabs((*p5)-(*p4)),(float)fabs((*p6)-(*p4)),(float)fabs((*p7)-(*p4)),(float)fabs((*p8)-(*p4)));

    return img_out;
}
