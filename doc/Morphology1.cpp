// MorphologyDlg.cpp : implementation file
//

#include "Morphology.h"

//形态学结构元素的复制
IplConvKernel* lhStructuringElementCopy(IplConvKernel* se)
{

	IplConvKernel* copy = cvCreateStructuringElementEx( se->nCols, se->nRows, 
			se->anchorX, se->anchorY, 0, NULL );

	copy->nShiftR = se->nShiftR;

	memcpy( copy->values, se->values, sizeof(int) * se->nRows * se->nCols );

	return copy;
}

//形态学结构元素的取反
IplConvKernel* lhStructuringElementNot(IplConvKernel* se)
{

	IplConvKernel* temp = cvCreateStructuringElementEx( se->nCols, se->nRows, 
			se->anchorX, se->anchorY, 0, NULL );

	temp->nShiftR = se->nShiftR;

	memcpy( temp->values, se->values, sizeof(int) * se->nRows * se->nCols );

	for(int i=0; i<temp->nRows * temp->nCols ; i++)
	{
		temp->values[i] = !temp->values[i] ;
	}

	return temp;
}
//形态学结构元素的基数
int lhStructuringElementCard(IplConvKernel* se)
{
	assert(se != NULL);
	int i, cnt = 0;

	for (i=0; i<se->nCols*se->nRows; i++)
	{
		cnt += se->values[i];
	}
	return cnt;

}

//形态学结构元素的映射
IplConvKernel* lhStructuringElementMap(IplConvKernel* se)
{
	CvMat *mat = cvCreateMat( se->nRows,  se->nCols, CV_32SC1);
	memcpy(mat->data.i, se->values, sizeof(int) * se->nRows * se->nCols );
	cvFlip(mat, NULL, -1);

	IplConvKernel* semap = cvCreateStructuringElementEx( se->nCols, se->nRows, 
			se->nCols-se->anchorX-1, se->nRows-se->anchorY-1, 0, NULL );

	semap->nShiftR = se->nShiftR;

	memcpy( semap->values, mat->data.i, sizeof(int) * se->nRows * se->nCols );

	cvReleaseMat(&mat);

	return semap;
}

//形态学线性结构元素的创建，常用于形态学方向分析
IplConvKernel* lhStructuringElementLine(unsigned int angle, unsigned int len)
{
	assert(len > 2);
	angle = angle%180;

	CvPoint p1 = {0};
	CvPoint p2 = {0};
	int width = cvRound(len*cos((float)angle*CV_PI/180.0));
	int height = cvRound(len*sin((float)angle*CV_PI/180.0));

	height = height >= 1 ? height : 1;

	if (width < 0)
	{
		width = width <= -1 ? width : -1;
		p1.x = 0;
		p1.y = 0;
		p2.x = -width -1;
		p2.y = height -1;
	}
	else
	{
		width = width >= 1 ? width : 1;
		p1.x = 0;
		p1.y = height -1;
		p2.x = width -1;
		p2.y = 0;
	}
	CvMat *temp = cvCreateMat(height, width, CV_32SC1);
	cvZero(temp);
	cvLine(temp, p1, p2, cvScalar(1,1,1), 1, 4, 0);

	IplConvKernel* se = cvCreateStructuringElementEx( width, height, 
			(width-1)/2, (height-1)/2,  CV_SHAPE_CUSTOM, temp->data.i );

	cvReleaseMat(&temp);
	return se;

}


//形态学开运算
void lhMorpOpen(const IplImage* src, IplImage* dst, IplConvKernel* se=NULL, int iterations=1)
{

    cvErode( src, dst, se, iterations );

	IplConvKernel* semap = lhStructuringElementMap(se);

    cvDilate( dst, dst, semap, iterations );

	cvReleaseStructuringElement(&semap);

}

//形态学闭运算
void lhMorpClose(const IplImage* src, IplImage* dst, IplConvKernel* se=NULL, int iterations=1)
{

    cvDilate( src, dst, se, iterations );

	IplConvKernel* semap = lhStructuringElementMap(se);

    cvErode( dst, dst, semap, iterations );

	cvReleaseStructuringElement(&semap);

}

//形态学基本梯度运算
void lhMorpGradient(const IplImage* src, IplImage* dst, IplConvKernel* se=NULL, int iterations=1)
{
	assert(src != NULL && dst != NULL && src != dst);
	IplImage*  temp = cvCloneImage(src);
	cvErode( src, temp, se, iterations );
    cvDilate( src, dst, se, iterations );
    cvSub( dst, temp, dst );
	cvReleaseImage(&temp);
}

//形态学内梯度运算
void lhMorpGradientIn(const IplImage* src, IplImage* dst, IplConvKernel* se=NULL, int iterations=1)
{
	assert(src != NULL && dst != NULL && src != dst);
	cvErode( src, dst, se, iterations );
    cvSub( src, dst, dst );
}

//形态学外梯度运算
void lhMorpGradientOut(const IplImage* src, IplImage* dst, IplConvKernel* se=NULL, int iterations=1)
{
	assert(src != NULL && dst != NULL && src != dst);
	cvDilate( src, dst, se, iterations );
    cvSub(dst, src, dst );
}

//形态学方向梯度
void lhMorpGradientDir(const IplImage* src, IplImage* dst, unsigned int angle, unsigned int len )
{
	assert(src != NULL && dst != NULL && src != dst);
	IplConvKernel* se = lhStructuringElementLine(angle, len);
	lhMorpGradient(src, dst, se);
	cvReleaseStructuringElement(&se);
}

//形态学白顶帽运算
void lhMorpWhiteTopHat(const IplImage* src, IplImage* dst, IplConvKernel* se=NULL, int iterations=1)
{
	assert(src != NULL && dst != NULL && src != dst);
	lhMorpOpen(src, dst, se, iterations );
    cvSub( src, dst, dst );
}


//形态学黑顶帽运算
void lhMorpBlackTopHat(const IplImage* src, IplImage* dst, IplConvKernel* se=NULL, int iterations=1)
{
	assert(src != NULL && dst != NULL && src != dst);
	lhMorpClose(src, dst, se, iterations );
    cvSub(dst, src, dst );

}

//形态学自补顶帽运算
void lhMorpQTopHat(const IplImage* src, IplImage* dst, IplConvKernel* se=NULL, int iterations=1)
{
	assert(src != NULL && dst != NULL && src != dst);
	IplImage*  temp = cvCloneImage(src);
    lhMorpClose( src, temp, se, iterations );
    lhMorpOpen( src, dst, se, iterations );
    cvSub(temp, dst, dst );
	cvReleaseImage(&temp);
}

//形态学对比度增强运算
void lhMorpEnhance(const IplImage* src, IplImage* dst, IplConvKernel* se=NULL, int iterations=1)
{
	assert(src != NULL && dst != NULL && src != dst);
	IplImage*  temp = cvCloneImage(src);
    lhMorpWhiteTopHat( src, temp, se, iterations );
    lhMorpBlackTopHat( src, dst, se, iterations );
	cvAdd(src, temp, temp);
    cvSub(temp, dst, dst );
	cvReleaseImage(&temp);
}

//形态学二值击中-击不中变换
void lhMorpHMTB(const IplImage* src, IplImage* dst, IplConvKernel* sefg, IplConvKernel* sebg =NULL)
{
	assert(src != NULL && dst != NULL && src != dst && sefg!= NULL && sefg!=sebg);

	if (sebg == NULL)
	{
		sebg = lhStructuringElementNot(sefg);

	}
	IplImage*  temp1 = cvCreateImage(cvGetSize(src), 8, 1);
	IplImage*  temp2 = cvCreateImage(cvGetSize(src), 8, 1);

	//P104 (5.2)
	cvErode( src, temp1, sefg);
	cvNot(src, temp2);
	cvErode( temp2, temp2, sebg);
	cvAnd(temp1, temp2, dst);


	cvReleaseImage(&temp1);
	cvReleaseImage(&temp2);

	cvReleaseStructuringElement(&sebg);

}


//形态学非约束击中-击不中变换 针对二值和灰度图像
void lhMorpHMTU(const IplImage* src, IplImage* dst, IplConvKernel* sefg, IplConvKernel* sebg =NULL)
{
	assert(src != NULL && dst != NULL && src != dst && sefg!= NULL && sefg!=sebg);

	if (sebg == NULL)
	{
		sebg = lhStructuringElementNot(sefg);

	}
	
	IplImage*  temp = cvCreateImage(cvGetSize(src), 8, 1);
	IplImage*  mask = cvCreateImage(cvGetSize(src), 8, 1);
	cvZero(mask);

	//P106 (5.4)
	cvErode( src, temp, sefg);
	cvDilate(src, dst, sebg);
	cvCmp(temp, dst, mask, CV_CMP_GT);

	cvSub(temp, dst, dst, mask);
	cvNot(mask, mask);
	cvSet(dst, cvScalar(0), mask);

	//cvCopy(dst, mask);
	//cvSet(dst, cvScalar(255), mask);
	cvReleaseImage(&mask);
	cvReleaseImage(&temp);

	cvReleaseStructuringElement(&sebg);
}

//形态学约束击中-击不中变换 针对二值和灰度图像
void lhMorpHMTC(const IplImage* src, IplImage* dst, IplConvKernel* sefg, IplConvKernel* sebg =NULL)
{
	assert(src != NULL && dst != NULL && src != dst && sefg!= NULL && sefg!=sebg);

	if (sebg == NULL)
	{
		sebg = lhStructuringElementNot(sefg);

	}
	
	IplImage*  temp1 = cvCreateImage(cvGetSize(src), 8, 1);
	IplImage*  temp2 = cvCreateImage(cvGetSize(src), 8, 1);
	IplImage*  temp3 = cvCreateImage(cvGetSize(src), 8, 1);
	IplImage*  temp4 = cvCreateImage(cvGetSize(src), 8, 1);

	IplImage*  mask1 = cvCreateImage(cvGetSize(src), 8, 1);
	IplImage*  mask2 = cvCreateImage(cvGetSize(src), 8, 1);
	IplImage*  mask3 = cvCreateImage(cvGetSize(src), 8, 1);
	IplImage*  mask4 = cvCreateImage(cvGetSize(src), 8, 1);

	cvZero(mask1);
	cvZero(mask2);
	cvZero(mask3);
	cvZero(mask4);

	cvZero(dst);

	//P107 (5.5)
	cvErode( src, temp1, sebg);
	cvDilate( src, temp2, sebg);
	cvErode( src, temp3, sefg);
	cvDilate( src, temp4, sefg);

	cvCmp(src, temp3, mask1, CV_CMP_EQ);
	cvCmp(temp2, src,  mask2, CV_CMP_LT);
	cvAnd(mask1, mask2, mask2);

	cvCmp(src, temp4, mask3 , CV_CMP_EQ);
	cvCmp(temp1, src, mask4 , CV_CMP_GT);
	cvAnd(mask3, mask4, mask4);

	cvSub(src, temp2, dst, mask2);
	cvSub(temp1, src, dst, mask4);




	cvReleaseImage(&mask1);
	cvReleaseImage(&mask2);
	cvReleaseImage(&mask3);
	cvReleaseImage(&mask4);

	cvReleaseImage(&temp1);
	cvReleaseImage(&temp2);
	cvReleaseImage(&temp3);
	cvReleaseImage(&temp4);

	cvReleaseStructuringElement(&sebg);

}

#define LH_MORP_TYPE_BINARY			0
#define LH_MORP_TYPE_UNCONSTRAIN	1
#define LH_MORP_TYPE_CONSTRAIN		2

//形态学约束击中-击不中变换
void lhMorpHMT(const IplImage* src, IplImage* dst, IplConvKernel* sefg, IplConvKernel* sebg =NULL, int type=LH_MORP_TYPE_BINARY)
{
	switch(type)
	{
	case LH_MORP_TYPE_BINARY:
		{
			lhMorpHMTB(src, dst, sefg, sebg);
			break;
		}

	case LH_MORP_TYPE_UNCONSTRAIN:
		{
			lhMorpHMTU(src, dst, sefg, sebg);
			break;
		}

	case LH_MORP_TYPE_CONSTRAIN:
		{
			lhMorpHMTC(src, dst, sefg, sebg);
			break;
		}
		
	default:
		break;

	}

}

//形态学击中-击不中开变换 
void lhMorpHMTOpen(const IplImage* src, IplImage* dst, IplConvKernel* sefg, IplConvKernel* sebg =NULL, int type=LH_MORP_TYPE_BINARY)
{
	assert(src != NULL && dst != NULL && src != dst && sefg!= NULL && sefg!=sebg);

	IplConvKernel* semap = lhStructuringElementMap(sefg);

	IplImage*  temp = cvCreateImage(cvGetSize(src), 8, 1);

	//P110 (5.8)
	lhMorpHMT(src, temp, sefg, sebg, type);
	cvDilate(temp, dst, semap);

	cvReleaseImage(&temp);
	cvReleaseStructuringElement(&semap);

}

//形态学细化运算
void lhMorpThin(const IplImage* src, IplImage* dst, IplConvKernel* sefg, IplConvKernel* sebg =NULL, int type=LH_MORP_TYPE_BINARY)
{

	assert(src != NULL && dst != NULL && src != dst && sefg!= NULL && sefg!=sebg);


	lhMorpHMT(src, dst, sefg, NULL, type);
	cvSub(src, dst, dst);

}

//形态学细化匹配运算
void lhMorpThinFit(const IplImage* src, IplImage* dst, IplConvKernel* sefg, IplConvKernel* sebg =NULL, int type=LH_MORP_TYPE_BINARY)
{

	assert(src != NULL && dst != NULL && src != dst && sefg!= NULL && sefg!=sebg);

	lhMorpHMTOpen(src, dst, sefg, NULL, type);
	cvSub(src, dst, dst);
}

//形态学粗化运算
void lhMorpThick(const IplImage* src, IplImage* dst, IplConvKernel* sefg, IplConvKernel* sebg =NULL, int type=LH_MORP_TYPE_BINARY)
{

	assert(src != NULL && dst != NULL && src != dst && sefg!= NULL && sefg!=sebg);


	lhMorpHMT(src, dst, sefg, NULL, type);
	cvAdd(src, dst, dst);

}

//形态学粗化不匹配运算
void lhMorpThickMiss(const IplImage* src, IplImage* dst, IplConvKernel* sefg, IplConvKernel* sebg =NULL, int type=LH_MORP_TYPE_BINARY)
{

	assert(src != NULL && dst != NULL && src != dst && sefg!= NULL && sefg!=sebg);

	lhMorpHMTOpen(src, dst, sefg, NULL, type);
	cvAdd(src, dst, dst);
}


//比较两个图像是否相同， 0 相同
int  lhImageCmp(const IplImage* img1, const IplImage* img2)
{
	assert(img1->width == img2->width && img1->height == img2->height && img1->imageSize == img2->imageSize );
	return memcmp(img1->imageData, img2->imageData, img1->imageSize);
}

//形态学测地膨胀和膨胀重建运算
void lhMorpRDilate(const IplImage* src, const IplImage* msk, IplImage* dst, IplConvKernel* se = NULL, int iterations=-1)
{

	assert(src != NULL && msk != NULL && dst != NULL && src != dst );

	if(iterations < 0)
	{
		//膨胀重建
		cvMin(src, msk, dst);
		cvDilate(dst, dst, se);
		cvMin(dst, msk, dst);

		IplImage*  temp1 = cvCreateImage(cvGetSize(src), 8, 1);
		//IplImage*  temp2 = cvCreateImage(cvGetSize(src), 8, 1);

		do
		{
			//record last result
			cvCopy(dst, temp1);
			cvDilate(dst, dst, se);
			cvMin(dst, msk, dst);
			//cvCmp(temp1, dst, temp2, CV_CMP_NE );

		}
		//while(cvSum(temp2).val[0] != 0);
		while(lhImageCmp(temp1, dst)!= 0);

		cvReleaseImage(&temp1);
		//cvReleaseImage(&temp2);

		return;	

	}
	else if (iterations == 0)
	{
		cvCopy(src, dst);
	}
	else
	{

		//普通测地膨胀 p136(6.1)
		cvMin(src, msk, dst);
		cvDilate(dst, dst, se);
		cvMin(dst, msk, dst);

		for(int i=1; i<iterations; i++)
		{
			cvDilate(dst, dst, se);
			cvMin(dst, msk, dst);

		}

	}
}

//形态学测地腐蚀和腐蚀重建运算
void lhMorpRErode(const IplImage* src,  const IplImage* msk, IplImage* dst, IplConvKernel* se = NULL, int iterations=-1)
{

	assert(src != NULL  && msk != NULL && dst != NULL && src != dst );

	if(iterations < 0)
	{
		//腐蚀重建
		cvMax(src, msk, dst);
		cvErode(dst, dst, se);
		cvMax(dst, msk, dst);

		IplImage*  temp1 = cvCreateImage(cvGetSize(src), 8, 1);
		//IplImage*  temp2 = cvCreateImage(cvGetSize(src), 8, 1);

		do
		{
			//record last result
			cvCopy(dst, temp1);
			cvErode(dst, dst, se);
			cvMax(dst, msk, dst);
			//cvCmp(temp1, dst, temp2, CV_CMP_NE);

		}
		//while(cvSum(temp2).val[0] != 0);
		while(lhImageCmp(temp1, dst)!= 0);

		cvReleaseImage(&temp1);
		//cvReleaseImage(&temp2);

		return;	

	}
	else if (iterations == 0)
	{
		cvCopy(src, dst);
	}
	else
	{
		//普通测地腐蚀 p137(6.2)
		cvMax(src, msk, dst);
		cvErode(dst, dst, se);
		cvMax(dst, msk, dst);

		for(int i=1; i<iterations; i++)
		{
			cvErode(dst, dst, se);
			cvMax(dst, msk, dst);
		}
	}
}

//形态学开重建
void lhMorpROpen(const IplImage* src, IplImage* dst, IplConvKernel* se = NULL, int iterations=1)
{
	assert(src != NULL  && dst != NULL && src != dst );

	//p155(6.16)
	IplImage*  temp = cvCreateImage(cvGetSize(src), 8, 1);
	cvErode(src, temp, se, iterations);
	lhMorpRDilate(temp, src, dst, se, -1);
	cvReleaseImage(&temp);

}

//形态学闭重建
void lhMorpRClose(const IplImage* src, IplImage* dst, IplConvKernel* se = NULL, int iterations=1)
{
	assert(src != NULL  && dst != NULL && src != dst );

	//p155(6.17)
	IplImage*  temp = cvCreateImage(cvGetSize(src), 8, 1);
	cvDilate(src, temp, se, iterations);
	lhMorpRErode(temp, src, dst, se, -1);
	cvReleaseImage(&temp);

}

//形态学白帽重建
void lhMorpRWTH(const IplImage* src, IplImage* dst, IplConvKernel* se = NULL, int iterations=1)
{
	assert(src != NULL  && dst != NULL && src != dst );
	//p156
	lhMorpROpen(src, dst, se, iterations);
	cvSub(src, dst, dst);
}

//形态学黑帽重建
void lhMorpRBTH(const IplImage* src, IplImage* dst, IplConvKernel* se = NULL, int iterations=1)
{
	assert(src != NULL  && dst != NULL && src != dst );
	//p156
	lhMorpRClose(src, dst, se, iterations);
	cvSub(dst, src, dst);
}


//形态学测地自对偶和自对偶重建运算
void lhMorpRSelfDual(const IplImage* src, const IplImage* msk, IplImage* dst, IplConvKernel* se = NULL, int iterations=-1)
{
	assert(src != NULL  && msk != NULL && dst != NULL && src != dst );

	//p140(6.7) p142 (6.10)
	IplImage*  temp1 = cvCreateImage(cvGetSize(src), 8, 1);
	IplImage*  temp2 = cvCreateImage(cvGetSize(src), 8, 1);

	cvZero(temp2);

	lhMorpRDilate(src, msk, temp1, se, iterations);

	lhMorpRErode(src, msk, dst, se, iterations);

	cvCmp(src, msk, temp2, CV_CMP_LE);

	cvCopy(temp1, dst, temp2);

	cvReleaseImage(&temp1);
	cvReleaseImage(&temp2);
}



//形态学区域极小值
void lhMorpRMin(const IplImage* src, IplImage* dst, IplConvKernel* se = NULL)
{
	assert(src != NULL  &&  dst != NULL && src != dst );

	//p149 (6.14)
	IplImage*  temp = cvCreateImage(cvGetSize(src), 8, 1);

	cvAddS(src, cvScalar(1), temp);
	
	lhMorpRErode(temp, src, dst, se);

	cvSub(dst, src, dst);

	cvReleaseImage(&temp);

}

//形态学区域极大值
void lhMorpRMax(const IplImage* src, IplImage* dst, IplConvKernel* se = NULL)
{
	assert(src != NULL  &&  dst != NULL && src != dst );

	//p149 (6.13)
	IplImage*  temp = cvCreateImage(cvGetSize(src), 8, 1);

	cvSubS(src, cvScalar(1), temp);
	
	lhMorpRDilate(temp, src, dst, se);

	cvSub(src, dst, dst);

	cvReleaseImage(&temp);

}

//形态学H极大值
void lhMorpHMax(const IplImage* src, IplImage* dst, unsigned char h, IplConvKernel* se = NULL)
{
	assert(src != NULL  &&  dst != NULL && src != dst );

	//p150
	IplImage*  temp = cvCreateImage(cvGetSize(src), 8, 1);

	cvSubS(src, cvScalar(h), temp);
	
	lhMorpRDilate(temp, src, dst, se);

	cvReleaseImage(&temp);

}


//形态学H极小值
void lhMorpHMin(const IplImage* src, IplImage* dst, unsigned char h, IplConvKernel* se = NULL)
{
	assert(src != NULL  &&  dst != NULL && src != dst );

	//p150
	IplImage*  temp = cvCreateImage(cvGetSize(src), 8, 1);

	cvAddS(src, cvScalar(h), temp);
	
	lhMorpRErode(temp, src, dst, se);

	cvReleaseImage(&temp);

}

//形态学H凹变换
void lhMorpHConcave(const IplImage* src, IplImage* dst, unsigned char h, IplConvKernel* se = NULL)
{
	assert(src != NULL  &&  dst != NULL && src != dst );

	//p150

	lhMorpHMin(src, dst, h, se);
	cvSub(dst, src, dst);

}

//形态学H凸变换
void lhMorpHConvex(const IplImage* src, IplImage* dst, unsigned char h, IplConvKernel* se = NULL)
{
	assert(src != NULL  &&  dst != NULL && src != dst );

	//p150

	lhMorpHMax(src, dst, h, se);
	cvSub(src, dst, dst);

}

//形态学扩展极大值
void lhMorpEMax(const IplImage* src, IplImage* dst, unsigned char h, IplConvKernel* se = NULL)
{
	assert(src != NULL  &&  dst != NULL && src != dst );

	//p150
	IplImage*  temp = cvCreateImage(cvGetSize(src), 8, 1);

	lhMorpHMax(src, temp, h, se);
	lhMorpRMax(temp, dst, se);

	cvReleaseImage(&temp);

}

//形态学扩展极小值
void lhMorpEMin(const IplImage* src, IplImage* dst, unsigned char h, IplConvKernel* se = NULL)
{
	assert(src != NULL  &&  dst != NULL && src != dst );

	//p150
	IplImage*  temp = cvCreateImage(cvGetSize(src), 8, 1);

	lhMorpHMin(src, temp, h, se);
	lhMorpRMin(temp, dst, se);

	cvReleaseImage(&temp);

}


//形态学等级滤波器（二值,默认SE为矩形3*3）
void lhMorpRankFilterB(const IplImage* src, IplImage* dst, IplConvKernel* se = NULL, unsigned int rank = 0)
{
	assert(src != NULL  &&  dst != NULL && src != dst );

	bool defaultse = false;
	int card;
	if (se == NULL)
	{
		card = 3*3;
		assert(rank >= 0 && rank <= card);
		se = cvCreateStructuringElementEx(3, 3, 1, 1, CV_SHAPE_RECT);
		defaultse = true;
	}
	else
	{
		card = lhStructuringElementCard(se);
		assert(rank >= 0 && rank <= card);
	}

	//default rank is median
	if (rank == 0)
		rank = card/2+1;

	IplConvKernel* semap =	lhStructuringElementMap(se);

	CvMat *semat = cvCreateMat(semap->nRows, semap->nCols, CV_32FC1);

	int i;
	for (i=0; i<semap->nRows*semap->nCols; i++)
	{
		semat->data.fl[i] = semap->values[i];
	}

	cvThreshold(src, dst, 0, 1, CV_THRESH_BINARY);
	IplImage *temp = cvCreateImage(cvGetSize(dst), 8, 1);

	cvFilter2D(dst, temp, semat, cvPoint(semap->anchorX, semap->anchorY));

	cvThreshold(temp, dst, card-rank, 255, CV_THRESH_BINARY);

	cvReleaseMat(&semat);
	cvReleaseStructuringElement(&semap);

	if (defaultse)
		cvReleaseStructuringElement(&se);	
	
	cvReleaseImage(&temp);

}

//形态学重建应用1：去除边界的连通区域
void lhMorpRemoveBoderObj(const IplImage* src, IplImage* dst)
{
	IplImage *temp = cvCloneImage(src);
	//double min, max;
	//cvMinMaxLoc(src, &min, &max);
	
	//标记图像
	cvRectangle(temp, cvPoint(3,3), cvPoint(temp->width-7, temp->height-7), CV_RGB(0,0,0), -1);

	//将原图像作为掩模图像
	lhMorpRDilate(temp, src, dst);

	cvReleaseImage(&temp);
	//cvSet((CvArr*)src, cvScalar(min), dst);
	//cvCopy(src, dst);
	cvSub(src, dst, dst);
}


//形态学重建应用2：空洞的填充
void lhMorpFillHole(const IplImage* src, IplImage* dst)
{
	IplImage *temp = cvCloneImage(src);
	double min, max;
	cvMinMaxLoc(src, &min, &max);
	//标记图像
	cvRectangle(temp, cvPoint(3,3), cvPoint(temp->width-7, temp->height-7), CV_RGB(max,max,max), -1);

	//将原图像作为掩模图像
	lhMorpRErode(temp, src, dst);

	cvReleaseImage(&temp);
	//cvSub(src, dst, dst);
}


//
//
///////////////////////////////////////////////////////////////////////////////
//// CMorphologyDlg dialog
//
//CMorphologyDlg::CMorphologyDlg(CWnd* pParent /*=NULL*/)
//	: CDialog(CMorphologyDlg::IDD, pParent)
//{
//	//{{AFX_DATA_INIT(CMorphologyDlg)
//	m_Shape = 0;
//	m_nBinTh = 180;
//	m_nMethod = 0;
//	//}}AFX_DATA_INIT
//	// Note that LoadIcon does not require a subsequent DestroyIcon in Win32
//	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
//
//	for (int i=0; i<PLOTNUM; i++)
//	{
//		m_pOrgImage[i] = NULL;
//	}
//
//	m_pKernel = NULL;
//
//	m_nCurImgFlag = 0;
//
//
//}
//
//void CMorphologyDlg::DoDataExchange(CDataExchange* pDX)
//{
//	CDialog::DoDataExchange(pDX);
//	//{{AFX_DATA_MAP(CMorphologyDlg)
//	DDX_Control(pDX, IDC_SPIN_ROW, m_SpinRow);
//	DDX_Control(pDX, IDC_SPIN_COL, m_SpinCol);
//	DDX_Control(pDX, IDC_SPIN_ANCHORY, m_SpinAnchorY);
//	DDX_Control(pDX, IDC_SPIN_ANCHORX, m_SpinAnchorX);
//	DDX_Control(pDX, IDC_EDIT_ROW, m_RowEdit);
//	DDX_Control(pDX, IDC_EDIT_COL, m_ColEdit);
//	DDX_Control(pDX, IDC_EDIT_ANCHORY, m_AnchorYEdit);
//	DDX_Control(pDX, IDC_EDIT_ANCHORX, m_AnchorXEdit);
//	DDX_Radio(pDX, IDC_RADIO_SHAPE1, m_Shape);
//	DDX_Text(pDX, IDC_EDIT_BINARYTH, m_nBinTh);
//	DDV_MinMaxUInt(pDX, m_nBinTh, 0, 255);
//	DDX_CBIndex(pDX, IDC_COMBO_MORP, m_nMethod);
//	//}}AFX_DATA_MAP
//}
//
//BEGIN_MESSAGE_MAP(CMorphologyDlg, CDialog)
//	//{{AFX_MSG_MAP(CMorphologyDlg)
//	ON_WM_PAINT()
//	ON_WM_QUERYDRAGICON()
//	ON_BN_CLICKED(IDC_BUTTON_OPENFILE, OnButtonOpenFile)
//	ON_BN_CLICKED(IDC_BUTTON_SAVEFILE, OnButtonSaveFile)
//	ON_WM_DESTROY()
//	ON_BN_CLICKED(IDC_BUTTON_RELOAD, OnButtonReload)
//	ON_BN_CLICKED(IDC_BUTTON_KERNEL_MODIFY, OnButtonKernelModify)
//	ON_BN_CLICKED(IDC_BUTTON_BINARY, OnButtonBinary)
//	ON_BN_CLICKED(IDC_BUTTON_MORPHOLOGY, OnButtonMorphology)
//	ON_BN_CLICKED(IDC_BUTTON_UNDO, OnButtonUndo)
//	ON_BN_CLICKED(IDC_BUTTON_NOT, OnButtonNot)
//	//}}AFX_MSG_MAP
//END_MESSAGE_MAP()
//
///////////////////////////////////////////////////////////////////////////////
//// CMorphologyDlg message handlers
//
//BOOL CMorphologyDlg::OnInitDialog()
//{
//	CDialog::OnInitDialog();
//
//	// Set the icon for this dialog.  The framework does this automatically
//	//  when the application's main window is not a dialog
//	SetIcon(m_hIcon, TRUE);			// Set big icon
//	SetIcon(m_hIcon, FALSE);		// Set small icon
//
//	
//	m_SpinRow.SetBuddy(&m_RowEdit);
//	m_SpinRow.SetRange(1.f, 300.0f);
//	m_SpinRow.SetPos(KERNEL_ROW);
//	m_SpinRow.SetDelta(1.0f);
//
//	m_SpinCol.SetBuddy(&m_ColEdit);
//	m_SpinCol.SetRange(1.f, 300.0f);
//	m_SpinCol.SetPos(KERNEL_COL);
//	m_SpinCol.SetDelta(1.0f);
//
//	m_SpinAnchorY.SetBuddy(&m_AnchorYEdit);
//	m_SpinAnchorY.SetRange(0.f, 300.0f);
//	m_SpinAnchorY.SetPos(ANCHOR_Y);
//	m_SpinAnchorY.SetDelta(1.0f);
//
//	m_SpinAnchorX.SetBuddy(&m_AnchorXEdit);
//	m_SpinAnchorX.SetRange(0.f, 300.0f);
//	m_SpinAnchorX.SetPos(ANCHOR_X);
//	m_SpinAnchorX.SetDelta(1.0f);
//
//
//	CStatic *pSt[PLOTNUM];
//	pSt[0] = (CStatic *)GetDlgItem(IDC_STATIC_IMAGE1);
//	pSt[1] = (CStatic *)GetDlgItem(IDC_STATIC_IMAGE2);
//	pSt[2] = (CStatic *)GetDlgItem(IDC_STATIC_IMAGE3);
//	pSt[3] = (CStatic *)GetDlgItem(IDC_STATIC_IMAGE4);	
//	
//	for (int i=0; i<PLOTNUM; i++)
//	{
//		RECT rect;
//		pSt[i]->GetWindowRect(&rect);
//		ScreenToClient(&rect);
//		m_FrameOrg[i].Create( rect.right-rect.left, rect.bottom-rect.top, 8);
//	}
//	
//
//
//	for (i=0; i<PLOTNUM; i++)
//	{
//		if (m_pOrgImage[i] == NULL)
//		{
//			m_pOrgImage[i] = cvCreateImage(cvSize(IMGWIDTH, IMGHEIGHT), 8, 1);
//			cvSetZero(m_pOrgImage[i] );
//		}
//
//	}
//	//cvLine(m_pOrgImage[3], cvPoint(100, 100), cvPoint(100, 100), CV_RGB(255, 255,255), 8);
//	//cvLine(m_pOrgImage[3], cvPoint(100, 101), cvPoint(100, 101), CV_RGB(255, 255,255));
//
//	OnButtonKernelModify();
//	
//
//	return TRUE;  // return TRUE  unless you set the focus to a control
//}
//
//// If you add a minimize button to your dialog, you will need the code below
////  to draw the icon.  For MFC applications using the document/view model,
////  this is automatically done for you by the framework.
//
//void CMorphologyDlg::OnPaint() 
//{
//	if (IsIconic())
//	{
//		CPaintDC dc(this); // device context for painting
//
//		SendMessage(WM_ICONERASEBKGND, (WPARAM) dc.GetSafeHdc(), 0);
//
//		// Center icon in client rectangle
//		int cxIcon = GetSystemMetrics(SM_CXICON);
//		int cyIcon = GetSystemMetrics(SM_CYICON);
//		CRect rect;
//		GetClientRect(&rect);
//		int x = (rect.Width() - cxIcon + 1) / 2;
//		int y = (rect.Height() - cyIcon + 1) / 2;
//
//		// Draw the icon
//		dc.DrawIcon(x, y, m_hIcon);
//	}
//	else
//	{
//		CPaintDC dc(this); // device context for painting
//		IplImage* ShownImage;
//		CStatic *pSt[PLOTNUM];
//		pSt[0] = (CStatic *)GetDlgItem(IDC_STATIC_IMAGE1);
//		pSt[1] = (CStatic *)GetDlgItem(IDC_STATIC_IMAGE2);
//		pSt[2] = (CStatic *)GetDlgItem(IDC_STATIC_IMAGE3);
//		pSt[3] = (CStatic *)GetDlgItem(IDC_STATIC_IMAGE4);
//
//
//		for (int i=0; i<PLOTNUM; i++)
//		{
//
//			ShownImage = m_FrameOrg[i].GetImage();
//			//ShownImage->origin = 1;
//			
//			cvResize(m_pOrgImage[i], ShownImage, CV_INTER_LINEAR );
//
//			//计算视频区域位置
//
//			RECT rect;
//			pSt[i]->GetWindowRect(&rect);
//			ScreenToClient(&rect);
//			m_FrameOrg[i].Show(dc.GetSafeHdc(), rect.left, rect.top, 
//				rect.right-1, rect.bottom-1, 0, 0 );
//
//		}
//
//		CDialog::OnPaint();
//	
//	}
//}
//
//// The system calls this to obtain the cursor to display while the user drags
////  the minimized window.
//HCURSOR CMorphologyDlg::OnQueryDragIcon()
//{
//	return (HCURSOR) m_hIcon;
//}
//
//void CMorphologyDlg::OnButtonOpenFile() 
//{
//	cvSetZero(m_pOrgImage[0]);
//	cvSetZero(m_pOrgImage[1]);
//
//	
//	//CFileDialog dlg(TRUE, NULL, GetCurrentPath()+"\\*.*", 
//	CFileDialog dlg(TRUE, NULL, "*.*", 
//		OFN_FILEMUSTEXIST|OFN_PATHMUSTEXIST|OFN_HIDEREADONLY,
//		"IMG files (*.bmp; *.jpg) |*.bmp;*.jpg||",NULL);
//	char title[]= {"载入图像"};
//	dlg.m_ofn.lpstrTitle= title;
//
//	if (dlg.DoModal() == IDOK) 
//	{
//		
//		IplImage* img = cvLoadImage(dlg.GetPathName(), 0);//CV_LOAD_IMAGE_GRAYSCALE);
//
//		m_strFilePath = dlg.GetPathName();
//		
//		SetDlgItemText(IDC_EDIT_FILEPATH, dlg.GetPathName());
//		
//		CString size;
//		size.Format("%d*%d", cvGetSize( img ).width, cvGetSize( img ).height );
//		SetDlgItemText(IDC_STATIC_IMGSIZE, size);
//
///*		if ( cvGetSize( img ).height != IMGHEIGHT  || cvGetSize( img ).width != IMGWIDTH  )
//		{
//
//			AfxMessageBox("图像大小不匹配：" + dlg.GetPathName());
//
//			cvReleaseImage(&img);
//
//			
//			return ;
//
//		}
//*/	
//		//cvResize(img, m_pOrgImage[0]);
//
//		ResizeImages(cvGetSize(img));
//		cvCopy(img, m_pOrgImage[0]);
//		cvCopy(img, m_pOrgImage[3]);
//
//		m_nCurImgFlag = 0;
//
//		//cvReleaseImage(&img);
//
//		//ProcessImage();
//	}
//	
//	UpdateImageRect();
//	
//}
//void CMorphologyDlg::ResizeImages(CvSize size)
//{
//	for (int i=0; i<PLOTNUM; i++)
//	{
//		if (m_pOrgImage[i] != NULL)
//		{
//			cvReleaseImage(&m_pOrgImage[i]);	
//			m_pOrgImage[i] = cvCreateImage(size, 8, 1);
//			cvSetZero(m_pOrgImage[i] );
//		}
//	}
//
//}
//
//
//void CMorphologyDlg::OnButtonSaveFile() 
//{
//	CFileDialog dlg(FALSE, ".bmp", GetCurrentPath()+"\\*.BMP", 
//		OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
//		"All Files (*.bmp)|*.bmp|",NULL);
//	char title[]= {"保存检测结果图像"};
//	dlg.m_ofn.lpstrTitle= title;
//
//	if (dlg.DoModal() == IDOK) 
//	{
//
//		cvSaveImage( dlg.GetPathName(), m_pOrgImage[3]);
//		SetDlgItemText(IDC_EDIT_FILEPATH2, dlg.GetPathName() );
//		
//	}
//	
//}
//
//CString CMorphologyDlg::GetCurrentPath()
//{
//
//	char buffer[MAX_PATH];
//	GetModuleFileName(NULL,buffer,MAX_PATH);
//	CString path =CString(buffer);
//	path = path.Left(path.ReverseFind('\\'));
//	return path;
//}
//
//void CMorphologyDlg::ProcessImage()
//{
//
//	int64 begin_count = cvGetTickCount();
//	
//	//cvZero(m_pOrgImage[1]);
//	//cvZero(m_pOrgImage[2]);
//	//cvZero(m_pOrgImage[3]);
//
//	cvCopyImage(m_pOrgImage[2], m_pOrgImage[1]);
//	cvCopyImage(m_pOrgImage[3], m_pOrgImage[2]);
//
//
//	IplImage *src = m_pOrgImage[2];
//	IplImage *dst = m_pOrgImage[3];
//
//
//
//
//	UpdateData(TRUE);
//	switch(m_nMethod)
//	{
//
//		//腐蚀
//		case 0:
//			{
//				cvErode(src, dst, m_pKernel);
//
//
//				break;
//			}
//
//		//膨胀
//		case 1:
//			{
//
//				cvDilate(src, dst, m_pKernel);
//
//				break;
//			}
//
//			
//		//开
//		case 2:
//			{
//				lhMorpOpen(src, dst, m_pKernel);
//
//
//				break;
//			}
//
//		//闭
//		case 3:
//			{
//
//				lhMorpClose(src, dst, m_pKernel);
//
//				break;
//			}
//
//		//形态梯度
//		case 4:
//			{
//
//				lhMorpGradient(src, dst, m_pKernel);
//
//				break;
//			}
//
//		//形态内梯度
//		case 5:
//			{
//
//				lhMorpGradientIn(src, dst, m_pKernel);
//
//				break;
//			}
//
//		//形态外梯度
//		case 6:
//			{
//
//				lhMorpGradientOut(src, dst, m_pKernel);
//
//				break;
//			}		
//			
//		//白帽
//		case 7:
//			{
//
//				lhMorpWhiteTopHat(src, dst, m_pKernel);
//
//				break;
//			}
//
//		//黑帽
//		case 8:
//			{
//
//				lhMorpBlackTopHat(src, dst, m_pKernel);
//
//				break;
//			}
//
//		//自补顶帽
//		case 9:
//			{
//
//				lhMorpQTopHat(src, dst, m_pKernel);
//
//				break;
//			}
//
//		//对比度增强
//		case 10:
//			{
//
//				lhMorpEnhance(src, dst, m_pKernel);
//
//				break;
//			}
//
//		//击中-击不中(二值)
//		case 11:
//			{
//
//				lhMorpHMT(src, dst, m_pKernel, LH_MORP_TYPE_BINARY);
//				break;
//			}
//
//
//		//击中-击不中(非约束)
//		case 12:
//			{
//
//				lhMorpHMT(src, dst, m_pKernel, NULL, LH_MORP_TYPE_UNCONSTRAIN);
//				//为了方便显示结果，二值化
//				cvThreshold(dst, dst, 0, 255, CV_THRESH_BINARY);
//
//				break;
//			}
//
//
//		//击中-击不中(约束)
//		case 13:
//			{
//
//				lhMorpHMT(src, dst, m_pKernel, NULL, LH_MORP_TYPE_CONSTRAIN);
//				//为了方便显示结果，二值化
//				cvThreshold(dst, dst, 0, 255, CV_THRESH_BINARY);
//
//				break;
//			}
//
//		//击中-击不中开(二值)
//		case 14:
//			{
//
//				lhMorpHMTOpen(src, dst, m_pKernel,  NULL, LH_MORP_TYPE_BINARY);
//
//				break;
//			}
//
//		//击中-击不中开(非约束)
//		case 15:
//			{
//
//				lhMorpHMTOpen(src, dst, m_pKernel,  NULL, LH_MORP_TYPE_UNCONSTRAIN);
//				//为了方便显示结果，二值化
//				cvThreshold(dst, dst, 0, 255, CV_THRESH_BINARY);
//				break;
//			}
//
//		//击中-击不中开(约束)
//		case 16:
//			{
//
//				lhMorpHMTOpen(src, dst, m_pKernel,  NULL, LH_MORP_TYPE_CONSTRAIN);
//				//为了方便显示结果，二值化
//				cvThreshold(dst, dst, 0, 255, CV_THRESH_BINARY);
//				break;
//			}
//
//		//细化(二值)
//		case 17:
//			{
//				lhMorpThin(src, dst, m_pKernel,  NULL, LH_MORP_TYPE_BINARY);
//
//				break;
//			}
//
//		//细化(非约束 灰度)
//		case 18:
//			{
//				lhMorpThin(src, dst, m_pKernel,  NULL, LH_MORP_TYPE_UNCONSTRAIN);
//
//				break;
//			}
//
//		//细化(约束 灰度)
//		case 19:
//			{
//				lhMorpThin(src, dst, m_pKernel,  NULL, LH_MORP_TYPE_CONSTRAIN);
//
//				break;
//			}
//
//		//细化匹配(二值)
//		case 20:
//			{
//				lhMorpThinFit(src, dst, m_pKernel,  NULL, LH_MORP_TYPE_BINARY);
//
//				break;
//			}
//
//		//细化匹配(非约束 灰度)
//		case 21:
//			{
//				lhMorpThinFit(src, dst, m_pKernel,  NULL, LH_MORP_TYPE_UNCONSTRAIN);
//
//				break;
//			}
//
//		//细化匹配(约束 灰度)
//		case 22:
//			{
//				lhMorpThinFit(src, dst, m_pKernel,  NULL, LH_MORP_TYPE_CONSTRAIN);
//
//				break;
//			}
//
//		//粗化(二值)
//		case 23:
//			{
//				lhMorpThick(src, dst, m_pKernel,  NULL, LH_MORP_TYPE_BINARY);
//
//				break;
//			}
//
//		//粗化(非约束 灰度)
//		case 24:
//			{
//				lhMorpThick(src, dst, m_pKernel,  NULL, LH_MORP_TYPE_UNCONSTRAIN);
//
//				break;
//			}
//
//		//粗化(约束 灰度)
//		case 25:
//			{
//				lhMorpThick(src, dst, m_pKernel,  NULL, LH_MORP_TYPE_CONSTRAIN);
//
//				break;
//			}
//
//		//粗化不匹配(二值)
//		case 26:
//			{
//				lhMorpThickMiss(src, dst, m_pKernel,  NULL, LH_MORP_TYPE_BINARY);
//
//				break;
//			}
//
//		//粗化不匹配(非约束 灰度)
//		case 27:
//			{
//				lhMorpThickMiss(src, dst, m_pKernel,  NULL, LH_MORP_TYPE_UNCONSTRAIN);
//
//				break;
//			}
//
//		//粗化不匹配(约束 灰度)
//		case 28:
//			{
//				lhMorpThickMiss(src, dst, m_pKernel,  NULL, LH_MORP_TYPE_CONSTRAIN);
//
//				break;
//			}
//
//
//		//开重建
//		case 29:
//			{
//
//				lhMorpROpen(src, dst, m_pKernel, 2);
//
//				break;
//			}
//
//		//闭重建
//		case 30:
//			{
//				lhMorpRClose(src, dst, m_pKernel, 2);
//
//				break;
//			}
//
//		//顶帽重建
//		case 31:
//			{
//
//				//白顶帽重建
//				//lhMorpRWTH(src, dst, m_pKernel, 2);
//				//黑顶帽重建
//				lhMorpRBTH(src, dst, m_pKernel, 2);
//				//为了方便显示结果，二值化
//				cvThreshold(dst, dst, 0, 255, CV_THRESH_BINARY);
//				break;
//			}
//
//		//区域极小值
//		case 32:
//			{
//				lhMorpRMin(src, dst, m_pKernel);
//				//为了方便显示结果，二值化
//				cvThreshold(dst, dst, 0, 255, CV_THRESH_BINARY);
//				break;
//			}
//
//		//区域极大值
//		case 33:
//			{
//
//				lhMorpRMax(src, dst, m_pKernel);
//				//为了方便显示结果，二值化
//				cvThreshold(dst, dst, 0, 255, CV_THRESH_BINARY);
//				break;
//			}
//
//		//H极小值
//		case 34:
//			{
//				lhMorpHMin(src, dst, 100, m_pKernel);
//
//				break;
//			}
//
//		//H极大值
//		case 35:
//			{
//
//				lhMorpHMax(src, dst, 100, m_pKernel);
//				break;
//			}
//
//		//H凹变换
//		case 36:
//			{
//				lhMorpHConcave(src, dst, 100, m_pKernel);
//				//为了方便显示结果，二值化
//				cvThreshold(dst, dst, 0, 255, CV_THRESH_BINARY);
//				break;
//			}
//
//		//H凸变换
//		case 37:
//			{
//				lhMorpHConvex(src, dst, 100, m_pKernel);
//				//为了方便显示结果，二值化
//				cvThreshold(dst, dst, 0, 255, CV_THRESH_BINARY);
//				break;
//			}
//
//		//扩展极小值
//		case 38:
//			{
//
//				lhMorpEMin(src, dst, 100, m_pKernel);
//				//为了方便显示结果，二值化
//				cvThreshold(dst, dst, 0, 255, CV_THRESH_BINARY);
//				break;
//			}
//
//		//扩展极大值
//		case 39:
//			{
//
//				lhMorpEMax(src, dst, 100, m_pKernel);
//				//为了方便显示结果，二值化
//				cvThreshold(dst, dst, 0, 255, CV_THRESH_BINARY);
//				break;
//			}
//		//等级滤波
//		case 40:
//			{
//
//				//中值滤波
//				//lhMorpRankFilterB(src, dst, m_pKernel);
//
//				//腐蚀
//				//lhMorpRankFilterB(src, dst, m_pKernel, 1);
//
//				//膨胀
//				lhMorpRankFilterB(src, dst, m_pKernel, lhStructuringElementCard(m_pKernel));
//				break;
//			}
//
//		//方向梯度 垂直
//		case 41:
//			{
//				//IplConvKernel* se = lhStructuringElementLine(45, 10);
//				//cvReleaseStructuringElement(&se);
//
//				lhMorpGradientDir(src, dst, 90, 3);
//
//				break;
//
//			}
//
//		//方向梯度 水平
//		case 42:
//			{
//				lhMorpGradientDir(src, dst, 0, 3);
//				/*
//				IplConvKernel* se = lhStructuringElementLine(0, 35);
//				lhMorpOpen(src, dst, se);
//				cvReleaseStructuringElement(&se);
//
//				IplImage *temp = cvCreateImage(cvGetSize(src), 8,1 );
//				se = lhStructuringElementLine(90, 35);
//				lhMorpOpen(src, temp, se);
//				cvReleaseStructuringElement(&se);				
//
//				cvSub(src, dst, dst);
//				cvSub(dst, temp, dst);
//				cvCop
//
//				lhMorpOpen(dst, dst);
//				*/
//				break;
//			}
//
//		//重建应用1：去除边界的连通区域
//		case 43:
//			{
//
//				lhMorpRemoveBoderObj(src, dst);
//				break;
//			}
//
//		//重建应用2：空洞的填充
//		case 44:
//			{
//
//				lhMorpFillHole(src, dst);
//				break;
//			}
//
//		default:
//			break;
//	}
//
//	UpdateImageRect();
//
//	SetDlgItemInt(IDC_STATIC_TIMECOST, (cvGetTickCount() - begin_count)/cvGetTickFrequency()/1000 );
// 
//}
//
//void CMorphologyDlg::OnDestroy() 
//{
//
//	for (int i=0; i<PLOTNUM; i++)
//	{
//		if (m_pOrgImage[i] != NULL)
//			cvReleaseImage(&m_pOrgImage[i]);	
//	}
//	if ( m_pKernel != NULL)
//	{
//		cvReleaseStructuringElement(&m_pKernel);
//		m_pKernel = NULL;
//	}
//	
//	CDialog::OnDestroy();
//	
//	
//}
//
//void CMorphologyDlg::UpdateImageRect()
//{
//	CStatic *pSt[PLOTNUM];
//	pSt[0] = (CStatic *)GetDlgItem(IDC_STATIC_IMAGE1);
//	pSt[1] = (CStatic *)GetDlgItem(IDC_STATIC_IMAGE2);
//	pSt[2] = (CStatic *)GetDlgItem(IDC_STATIC_IMAGE3);
//	pSt[3] = (CStatic *)GetDlgItem(IDC_STATIC_IMAGE4);
//
//	for (int i=0; i<PLOTNUM; i++)
//	{
//		RECT rect;
//		pSt[i]->GetWindowRect(&rect);
//		ScreenToClient(&rect);
//		InvalidateRect(&rect, FALSE);
//	}
//
//}
//
//void CMorphologyDlg::OnButtonReload() 
//{
//	if (m_strFilePath.IsEmpty())
//		return;
//
//	IplImage* img = cvLoadImage(m_strFilePath, 0);//CV_LOAD_IMAGE_GRAYSCALE);
//
//	ResizeImages(cvGetSize(img));
//	cvCopy(img, m_pOrgImage[0]);
//	cvCopy(img, m_pOrgImage[3]);
//
//	m_nCurImgFlag = 0;
//
//	UpdateImageRect();	
//}
//
//void CMorphologyDlg::OnButtonKernelModify() 
//{
//	UpdateData(TRUE);
//
//	if ( m_pKernel != NULL)
//	{
//		cvReleaseStructuringElement(&m_pKernel);
//		m_pKernel = NULL;
//	}
//
//	//防止anchor值溢出
//	if(m_SpinAnchorX.GetPos() >= m_SpinCol.GetPos())
//	{
//		m_SpinAnchorX.SetPos((int)m_SpinCol.GetPos()/2);
//
//	}
//
//	if(m_SpinAnchorY.GetPos() >= m_SpinRow.GetPos())
//	{
//		m_SpinAnchorY.SetPos((int)m_SpinRow.GetPos()/2);
//
//	}
//
//	if(m_Shape != 3)
//	{
//
//		m_pKernel = cvCreateStructuringElementEx( m_SpinCol.GetPos(), m_SpinRow.GetPos(), 
//			m_SpinAnchorX.GetPos(), m_SpinAnchorY.GetPos(), m_Shape );
//	}
//	else
//	{
//		//自定义结构元素SE
//		int col = 9, row = 4;
//		int anchorx = 1, anchory = 0;
//		int kernel[] = {0, 0, 0, 1, 1, 1, 0, 0, 0 };
//		//int kernel[] = {0, 1};
//
//		m_pKernel = cvCreateStructuringElementEx( col, row, anchorx, anchory, CV_SHAPE_CUSTOM, kernel );
//
//		//更新界面参数显示
//		m_SpinCol.SetPos(col);
//		m_SpinRow.SetPos(row);
//		m_SpinAnchorX.SetPos(anchorx);
//		m_SpinAnchorY.SetPos(anchory);
//		//UpdateData(FALSE);
//	}
//
//	//int temp[256];
//	//memcpy(temp, m_pKernel->values, m_SpinCol.GetPos()*m_SpinRow.GetPos()*sizeof(int));
//
//}
//
//void CMorphologyDlg::OnButtonBinary() 
//{
//	UpdateData(TRUE);
//
//	cvCopyImage(m_pOrgImage[2], m_pOrgImage[1]);
//	cvCopyImage(m_pOrgImage[3], m_pOrgImage[2]);
//
//
//	IplImage *src = m_pOrgImage[2];
//	IplImage *dst = m_pOrgImage[3];
//	
//	cvThreshold(src, dst, m_nBinTh, 255, CV_THRESH_BINARY);
//	UpdateImageRect();
//	
//}
//
//void CMorphologyDlg::OnButtonMorphology() 
//{
//	IplImage*  temp1 = cvCreateImage(cvSize(200, 200), IPL_DEPTH_8U, 1);
//	IplImage*  temp2 = cvCreateImage(cvSize(200, 200), IPL_DEPTH_8U, 1);
//
//	cvSet(temp1, cvScalar(1));
//	cvSet(temp2, cvScalar(255));
//	cvSub(temp1, temp2, temp2);
//
//	cvReleaseImage(&temp1);
//	cvReleaseImage(&temp2);
//
//	ProcessImage();
//	
//}
//
//void CMorphologyDlg::OnButtonUndo() 
//{
//	cvCopyImage(m_pOrgImage[2], m_pOrgImage[3]);
//	cvCopyImage(m_pOrgImage[1], m_pOrgImage[2]);
//	UpdateImageRect();
//	
//}
//
//void CMorphologyDlg::OnButtonNot() 
//{
//	cvCopyImage(m_pOrgImage[2], m_pOrgImage[1]);
//	cvCopyImage(m_pOrgImage[3], m_pOrgImage[2]);
//
//
//	IplImage *src = m_pOrgImage[2];
//	IplImage *dst = m_pOrgImage[3];
//	
//	cvNot(src, dst);
//	UpdateImageRect();
//	
//}
