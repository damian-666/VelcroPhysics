using System;
using System.Collections.Generic;
using System.Text;

namespace IWAVE
{
    class IWave
    {

        /*
         #include <cmath>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/gl.h> // OpenGL itself.
#include <GL/glu.h> // GLU support library.
#include <GL/glut.h> // GLUT support library.
#endif
#include <iostream>
using namespace std;
int iwidth, iheight, size;
float *display_map;
float *obstruction;
float *source;
float *height;
float *previous_height;
float *vertical_derivative;
float scaling_factor;
float kernel[13][13];
int paint_mode;
enum{ PAINT_OBSTRUCTION, PAINT_SOURCE };
bool regenerate_data;
bool toggle_animation_on_off;
float dt, alpha, gravity;
float obstruction_brush[3][3];
float source_brush[3][3];
int xmouse_prev, ymouse_prev;
//--------------------------------------------------------
//
// Initialization routines
//
//
// Initialize all of the fields to zero
void Initialize( float *data, int size, float value )
{
    for(int i=0;i<size;i++ ) { data[i] = value; }
}
// Compute the elements of the convolution kernel
void InitializeKernel()
{
    double dk = 0.01;
    double sigma = 1.0;
    double norm = 0;
    for(double k=0;k<10;k+=dk)
    {
        norm += k*k*exp(-sigma*k*k);
    }
    for( int i=-6;i<=6;i++ )
    {
        for( int j=-6;j<=6;j++ )
        {
            double r = sqrt( (float)(i*i + j*j) );
            double kern = 0;
            for( double k=0;k<10;k+=dk)
            {
                kern += k*k*exp(-sigma*k*k)*j0(r*k);
            }
            kernel[i+6][j+6] = kern / norm;
        }
    }
}
void InitializeBrushes()
{
    obstruction_brush[1][1] = 0.0;
    obstruction_brush[1][0] = 0.5;
    obstruction_brush[0][1] = 0.5;
    obstruction_brush[2][1] = 0.5;
    obstruction_brush[1][2] = 0.5;
    obstruction_brush[0][2] = 0.75;
    obstruction_brush[2][0] = 0.75;
    obstruction_brush[0][0] = 0.75;
    obstruction_brush[2][2] = 0.75;
    source_brush[1][1] = 1.0;
    source_brush[1][0] = 0.5;
    source_brush[0][1] = 0.5;
    source_brush[2][1] = 0.5;
    source_brush[1][2] = 0.5;
    source_brush[0][2] = 0.25;
    source_brush[2][0] = 0.25;
    source_brush[0][0] = 0.25;
    source_brush[2][2] = 0.25;
}
void ClearObstruction()
{
    for(int i=0;i<size;i++ ){ obstruction[i] = 1.0; }
}
void ClearWaves()
{
    for(int i=0;i<size;i++ )
    {
        height[i] = 0.0;
        previous_height[i] = 0.0;
        vertical_derivative[i] = 0.0;
    }
}
//----------------------------------------------------

//----------------------------------------------------
//
// These two routines,
//
// ComputeVerticalDerivative()
// Propagate()
//
// are the heart of the iWave algorithm.
//
// In Propagate(), we have not bothered to handle the
// boundary conditions. This makes the outermost
// 6 pixels all the way around act like a hard wall.
//
void ComputeVerticalDerivative()
{
    // first step: the interior
    for(int ix=6;ix<iwidth-6;ix++)
    {
        for(int iy=6;iy<iheight-6;iy++)
        {
            int index = ix + iwidth*iy;
            float vd = 0;
            for(int iix=-6;iix<=6;iix++)
            {
                for(int iiy=-6;iiy<=6;iiy++)
                {
                    int iindex = ix+iix + iwidth*(iy+iiy);
                    vd += kernel[iix+6][iiy+6] * height[iindex];
                }
            }
            vertical_derivative[index] = vd;
        }
    }
}
void Propagate()
{
    // apply obstruction
    for( int i=0;i<size;i++ ) { height[i] *= obstruction[i]; }
    // compute vertical derivative
    ComputeVerticalDerivative();
    // advance surface
    float adt = alpha*dt;
    float adt2 = 1.0/(1.0+adt);
    for( int i=0;i<size;i++ )
    {
        float temp = height[i];
        height[i] = height[i]*(2.0-adt)-previous_height[i]-gravity*vertical_derivative[i];
        height[i] *= adt2;
        height[i] += source[i];
        height[i] *= obstruction[i];
        previous_height[i] = temp;
        // reset source each step
        source[i] = 0;
    }
}
//------------------------------------------
//
// Painting and display code
//
void resetScaleFactor( float amount )
{
    scaling_factor *= amount;
}
void DabSomePaint( int x, int y )
{
    int xstart = x - 1;
    int ystart = y - 1;
    if( xstart < 0 ){ xstart = 0; }
    if( ystart < 0 ){ ystart = 0; }
    int xend = x + 1;
    int yend = y + 1;
    if( xend >= iwidth ){ xend = iwidth-1; }
    if( yend >= iheight ){ yend = iheight-1; }
    if( paint_mode == PAINT_OBSTRUCTION )
    {
        for(int ix=xstart;ix <= xend; ix++)
        {
            for( int iy=ystart;iy<=yend; iy++)
            {
                int index = ix + iwidth*(iheight-iy-1);
                obstruction[index] *= obstruction_brush[ix-xstart][iy-ystart];
            }
        }
    }
    else if( paint_mode == PAINT_SOURCE )
    {
        for(int ix=xstart;ix <= xend; ix++)
        {
            for( int iy=ystart;iy<=yend; iy++)
            {
                int index = ix + iwidth*(iheight-iy-1);
                source[index] += source_brush[ix-xstart][iy-ystart];
            }
        }
    }
    return;
}
//----------------------------------------------------
//
// GL and GLUT callbacks
//
//----------------------------------------------------
void cbDisplay( void )
{
    glClear(GL_COLOR_BUFFER_BIT );
    glDrawPixels( iwidth, iheight, GL_LUMINANCE, GL_FLOAT, display_map );
    glutSwapBuffers();
}
// animate and display new result
void cbIdle()
{
    if( toggle_animation_on_off ) { Propagate(); }
    ConvertToDisplay();
    cbDisplay();
}
void cbOnKeyboard( unsigned char key, int x, int y )
{
    switch (key)
    {
    case ’-’: case ’_’:
        resetScaleFactor( 1.0/0.9 );
        regenerate_data = true;
        break;
    case ’+’: case ’=’:
        resetScaleFactor( 0.9 );
        regenerate_data = true;
        break;
    case ’ ’:
        toggle_animation_on_off = !toggle_animation_on_off;
    case ’o’:
        paint_mode = PAINT_OBSTRUCTION;
        break;
    case ’s’:
        paint_mode = PAINT_SOURCE;
        break;
    case ’b’:
        ClearWaves();
        ClearObstruction();
        Initialize( source, size, 0.0 );
        break;
    default:
        break;
    }
}
void cbMouseDown( int button, int state, int x, int y )
{
    if( button != GLUT_LEFT_BUTTON ) { return; }
    if( state != GLUT_DOWN ) { return; }
    xmouse_prev = x;
    ymouse_prev = y;
    DabSomePaint( x, y );
}
void cbMouseMove( int x, int y )
{
    xmouse_prev = x;
    ymouse_prev = y;
    DabSomePaint( x, y );
}
//---------------------------------------------------
int main(int argc, char** argv)
{
    // initialize a few variables
    iwidth = iheight = 200;
    size = iwidth*iheight;
    dt = 0.03;
    alpha = 0.3;
    gravity = 9.8 * dt * dt;
    scaling_factor = 1.0;
    toggle_animation_on_off = true;
    // allocate space for fields and initialize them
    height = new float[size];
    previous_height = new float[size];
    vertical_derivative = new float[size];
    obstruction = new float[size];
    source = new float[size];
    display_map = new float[size];
    ClearWaves();
    ClearObstruction();
    ConvertToDisplay();
    Initialize( source, size, 0 );
    InitializeBrushes();
    // build the convolution kernel
    InitializeKernel();
    // GLUT routines
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowSize( iwidth, iheight );
    // Open a window
    char title[] = "iWave Demo";
    int Window_ID = glutCreateWindow( title );
    glClearColor( 1,1,1,1 );
    glutDisplayFunc(&cbDisplay);
    glutIdleFunc(&cbIdle);
    glutKeyboardFunc(&cbOnKeyboard);
    glutMouseFunc( &cbMouseDown );
    glutMotionFunc( &cbMouseMove );
    glutMainLoop();
    return 1;
}; 
        */

    //    https://www.atnf.csiro.au/computing/software/gipsy/sub/bessel.c

        static double bessj0(double x)
        /*------------------------------------------------------------*/
        /* PURPOSE: Evaluate Bessel function of first kind and order  */
        /*          0 at input x                                      */
        /*------------------------------------------------------------*/
        {
            double ax, z;
            double xx, y, ans, ans1, ans2;

            if ((ax = Math.Abs(x)) < 8.0)
            {
                y = x * x;
                ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7
                   + y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
                ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718
                   + y * (59272.64853 + y * (267.8532712 + y * 1.0))));
                ans = ans1 / ans2;
            }
            else
            {
                z = 8.0 / ax;
                y = z * z;
                xx = ax - 0.785398164;
                ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
                   + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
                ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
                   + y * (-0.6911147651e-5 + y * (0.7621095161e-6
                   - y * 0.934935152e-7)));
                ans = Math.Sqrt(0.636619772 / ax) * (Math.Cos(xx) * ans1 - z * Math.Sin(xx) * ans2);
            }
            return ans;
        }





        int iwidth;
        int iheight;
        int size;




        private float[] display_map;
        private float[] obstruction;
        private float[] source;
          float[] height;
          float[] previous_height;
          float[] vertical_derivative;
        private float scaling_factor;
        float[,] kernel = new float[13, 13];


        enum ePaintType { PAINT_OBSTRUCTION, PAINT_SOURCE };
        ePaintType paint_mode;
        bool regenerate_data;
        public bool toggle_animation_on_off;
        float dt, alpha, gravity;
        float[,] obstruction_brush = new float[3, 3];
        float[,] source_brush = new float[3, 3];
        int xmouse_prev, ymouse_prev;

        //--------------------------------------------------------
        //
        // Initialization routines
        //
        //
        // Initialize all of the fields to zero


/* unsafe
int iwidth, iheight, size;
float *display_map;
float *obstruction;
float *source;
float *height;
float *previous_height;
float *vertical_derivative;
float scaling_factor;
float[][] kernel = new float[13][13];
int paint_mode;
enum ePaintType{ PAINT_OBSTRUCTION, PAINT_SOURCE };
bool regenerate_data;
bool toggle_animation_on_off;
float dt, alpha, gravity;
float[][] obstruction_brush= new float [3][3];
float[][]  source_brush = new float[3][3];
int xmouse_prev, ymouse_prev;
*/
//--------------------------------------------------------
//
// Initialization routines
//
//
// Initialize all of the fields to  UNzero





#if UNSAFE
        void Initialize( float *data, int size, float value )
{
    for(int i=0;i<size;i++ ) { data[i] = value; }
}
// Compute the elements of the convolution kernel



#endif
void InitializeKernel()
{
    double dk = 0.01;
    double sigma = 1.0;
    double norm = 0;
    for(double k=0;k<10;k+=dk)
    {
        norm += k*k* Math.Exp(-sigma*k*k);
    }
    for( int i=-6;i<=6;i++ )
    {
        for( int j=-6;j<=6;j++ )
        {
            double r = Math.Sqrt( (float)(i*i + j*j) );
            double kern = 0;
            for( double k=0;k<10;k+=dk)
            {
                kern += k*k*Math.Exp(-sigma*k*k)* bessj0(r*k);
            }
            kernel[i+6,j+6] = (float)(kern / norm);
        }
    }
}
void InitializeBrushes()
{
    obstruction_brush[1,1] = 0.0f;
    obstruction_brush[1,0] = 0.5f;
    obstruction_brush[0,1] = 0.5f;
    obstruction_brush[2,1] = 0.5f;
    obstruction_brush[1,2] = 0.5f;
    obstruction_brush[0,2] = 0.75f;
    obstruction_brush[2,0] = 0.75f;
    obstruction_brush[0,0] = 0.75f;
    obstruction_brush[2,2] = 0.75f;
    source_brush[1,1] = 1.0f;
    source_brush[1,0] = 0.5f;
    source_brush[0,1] = 0.5f;
    source_brush[2,1] = 0.5f;
    source_brush[1,2] = 0.5f;
    source_brush[0,2] = 0.25f;
    source_brush[2,0] = 0.25f;
    source_brush[0,0] = 0.25f;
    source_brush[2,2] = 0.25f;
}
void ClearObstruction()
{
    for(int i=0;i<size;i++ ){ obstruction[i] = 1.0f; }
}
void ClearWaves()
{
    for(int i=0;i<size;i++ )
    {
        height[i] = 0.0f;
        previous_height[i] = 0.0f;
        vertical_derivative[i] = 0.0f;
    }
}
//----------------------------------------------------
public void ConvertToDisplay()
{
    for(int i=0;i<size;i++ )
    {
        display_map[i] = 0.5f*( height[i]/scaling_factor + 1.0f )*obstruction[i];
    }
}
//----------------------------------------------------
//
// These two routines,
//
// ComputeVerticalDerivative()
// Propagate()
//
// are the heart of the iWave algorithm.
//
// In Propagate(), we have not bothered to handle the
// boundary conditions. This makes the outermost
// 6 pixels all the way around act like a hard wall.
//
void ComputeVerticalDerivative()
{
    // first step: the interior
    for(int ix=6;ix<iwidth-6;ix++)
    {
        for(int iy=6;iy<iheight-6;iy++)
        {
            int index = ix + iwidth*iy;
            float vd = 0;
            for(int iix=-6;iix<=6;iix++)
            {
                for(int iiy=-6;iiy<=6;iiy++)
                {
                    int iindex = ix+iix + iwidth*(iy+iiy);
                    vd += kernel[iix+6,iiy+6] * height[iindex];
                }
            }
            vertical_derivative[index] = vd;
        }
    }
}


 
        public void Propagate()
{
    // apply obstruction
    for( int i=0;i<size;i++ ) { height[i] *= obstruction[i]; }
    // compute vertical derivative
    ComputeVerticalDerivative();
    // advance surface
    float adt = alpha*dt;
    float adt2 = 1.0f/(1.0f+adt);
    for( int i=0;i<size;i++ )
    {
        float temp = height[i];
        height[i] = height[i]*(2.0f-adt)-previous_height[i]-gravity*vertical_derivative[i];
        height[i] *= adt2;
        height[i] += source[i];
        height[i] *= obstruction[i];
        previous_height[i] = temp;
        // reset source each step
        source[i] = 0;
    }
}
//------------------------------------------
//
// Painting and display code
//
void resetScaleFactor( float amount )
{
    scaling_factor *= amount;
}
public void DabSomePaint( int x, int y )
{
    int xstart = x - 1;
    int ystart = y - 1;
    if( xstart < 0 ){ xstart = 0; }
    if( ystart < 0 ){ ystart = 0; }
    int xend = x + 1;
    int yend = y + 1;
    if( xend >= iwidth ){ xend = iwidth-1; }
    if( yend >= iheight ){ yend = iheight-1; }
    if( paint_mode == ePaintType.PAINT_OBSTRUCTION)
    {
        for(int ix=xstart;ix <= xend; ix++)
        {
            for( int iy=ystart;iy<=yend; iy++)
            {
                int index = ix + iwidth*(iheight-iy-1);
                        obstruction[index] *= obstruction_brush[ix - xstart,iy - ystart];
            }
        }
    }
    else if( paint_mode ==  ePaintType.PAINT_SOURCE )
    {
        for(int ix=xstart;ix <= xend; ix++)
        {
            for( int iy=ystart;iy<=yend; iy++)
            {
                int index = ix + iwidth*(iheight-iy-1);
                source[index] += source_brush[ix-xstart,iy-ystart];
            }
        }
    }
    return;
}


#if PORT
//----------------------------------------------------
//
// GL and GLUT callbacks
//
//----------------------------------------------------
void cbDisplay( void )
{
    glClear(GL_COLOR_BUFFER_BIT );
    glDrawPixels( iwidth, iheight, GL_LUMINANCE, GL_FLOAT, display_map );
    glutSwapBuffers();
}
// animate and display new result
void cbIdle()
{
    if( toggle_animation_on_off ) { Propagate(); }
    ConvertToDisplay();
    cbDisplay();
}
#endif

#if PORT
void cbOnKeyboard( unsigned char key, int x, int y )
{
    switch (key)
    {
    case ’-’: case ’_’:
        resetScaleFactor( 1.0/0.9 );
        regenerate_data = true;
        break;
    case ’+’: case ’=’:
        resetScaleFactor( 0.9 );
        regenerate_data = true;
        break;
    case ’ ’:
        toggle_animation_on_off = !toggle_animation_on_off;
    case ’o’:
        paint_mode = PAINT_OBSTRUCTION;
        break;
    case ’s’:
        paint_mode = PAINT_SOURCE;
        break;
    case ’b’:
        ClearWaves();
        ClearObstruction();
        Initialize( source, size, 0.0 );
        break;
    default:
        break;
    }
}
void cbMouseDown( int button, int state, int x, int y )
{
    if( button != GLUT_LEFT_BUTTON ) { return; }
    if( state != GLUT_DOWN ) { return; }
    xmouse_prev = x;
    ymouse_prev = y;
    DabSomePaint( x, y );
}
void cbMouseMove( int x, int y )
{
    xmouse_prev = x;
    ymouse_prev = y;
    DabSomePaint( x, y );
}
#endif
//int main(int argc, char** argv)



       public void StartiWave()
        {
        // initialize a few variables
         iwidth = 200;
           
             iheight = 200;

    size = iwidth* iheight;
        dt = 0.03f;
    alpha = 0.3f;
    gravity = 9.8f * dt* dt;
        scaling_factor = 1.0f;
    toggle_animation_on_off = true;
    // allocate space for fields and initialize them
    height = new float[size];
    previous_height = new float[size];
    vertical_derivative = new float[size];
    obstruction = new float[size];
    source = new float[size];
    display_map = new float[size];
    ClearWaves();
        ClearObstruction();
        ConvertToDisplay();
        Initialize(source, size, 0 );
        InitializeBrushes();
        // build the convolution kernel
        InitializeKernel();
        // GLUT routines
    //    glutInit(&argc, argv);
   ////     glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
     //   glutInitWindowSize(iwidth, iheight );
        // Open a window
   //     char title[] = "iWave Demo";
    //    int Window_ID = glutCreateWindow(title);
    //    glClearColor( 1,1,1,1 );
   //     glutDisplayFunc(&cbDisplay);
   //     glutIdleFunc(&cbIdle);
   //     glutKeyboardFunc(&cbOnKeyboard);
    //    glutMouseFunc( &cbMouseDown );
    //    glutMotionFunc( &cbMouseMove );
   //     glutMainLoop();

        }

         void Initialize(float[] data, int size, float value)
        {
            for (int i = 0; i < size; i++) { data[i] = value; }
        }
        // Compute the elements of the convolution kernel
      
        //----------------------------------------------------
        //
        // GL and GLUT callbacks
        //
        //----------------------------------------------------
  /*      void cbDisplay(void )
        {
            glClear(GL_COLOR_BUFFER_BIT);
            glDrawPixels(iwidth, iheight, GL_LUMINANCE, GL_FLOAT, display_map);
            glutSwapBuffers();
        }
        // animate and display new result
        void cbIdle()
        {
            if (toggle_animation_on_off) { Propagate(); }
            ConvertToDisplay();
            cbDisplay();
        }


        void cbOnKeyboard(unsigned char key, int x, int y)
        {
            switch (key)
            {
                case ’-’: case ’_’:
        resetScaleFactor(1.0 / 0.9);
                    regenerate_data = true;
                    break;
                case ’+’: case ’=’:
        resetScaleFactor(0.9);
                    regenerate_data = true;
                    break;
                case ’ ’:
        toggle_animation_on_off = !toggle_animation_on_off;
                case ’o’:
        paint_mode = PAINT_OBSTRUCTION;
                    break;
                case ’s’:
        paint_mode = PAINT_SOURCE;
                    break;
                case ’b’:
        ClearWaves();
                    ClearObstruction();
                    Initialize(source, size, 0.0);
                    break;
                default:
                    break;
            }
        }
        void cbMouseDown(int button, int state, int x, int y)
        {
            if (button != GLUT_LEFT_BUTTON) { return; }
            if (state != GLUT_DOWN) { return; }
            xmouse_prev = x;
            ymouse_prev = y;
            DabSomePaint(x, y);
        }
        void cbMouseMove(int x, int y)
        {
            xmouse_prev = x;
            ymouse_prev = y;
            DabSomePaint(x, y);
        }

        */
        //---------------------------------------------------
         int DisplayWave()
        {
            // initialize a few variables
            iwidth = iheight = 200;
            size = iwidth * iheight;
            dt = 0.03f;
            alpha = 0.3f;
            gravity = 9.8f * dt * dt;
            scaling_factor = 1.0f;
            toggle_animation_on_off = true;
            // allocate space for fields and initialize them

            height = new float[size];
            previous_height = new float[size];
            vertical_derivative = new float[size];
            obstruction = new float[size];
            source = new float[size];
            display_map = new float[size];


            ClearWaves();
            ClearObstruction();
            ConvertToDisplay();
            Initialize(source, size, 0);
            InitializeBrushes();
            // build the convolution kernel
            InitializeKernel();
            // GLUT routines
       
      //      glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
           // glutInitWindowSize(iwidth, iheight);

            // Open a window
         // string title= "iWave Demo";
         //   int Window_ID = glutCreateWindow(title);
          //  glClearColor(1, 1, 1, 1);
         //   glutDisplayFunc(&cbDisplay);
       //     glutIdleFunc(&cbIdle);
       //     glutKeyboardFunc(&cbOnKeyboard);
       //     glutMouseFunc(&cbMouseDown);
     //       glutMotionFunc(&cbMouseMove);
    //        glutMainLoop();
            return 1;
        }
    }
}



/*
 * 
 * some more cluies fomr a supposed only tranfer of the book game-programming-gems-4
3.6 


Interactive Water Surfaces 

Jerry Tessendorf, 

Rhythm & Hues Studios 

jerryt@rhythm.com 


R ealistic computer-generated ocean surfaces have been used routinely in film 
features since 1996, in such titles as Waterworld, Titanic, Fifth Element, The 
Perfect Storm, X2 XMen United, Finding Nemo, and many more. For the most part, 
the algorithms underlying these productions apply Fast Fourier Transforms (FFT) 
to carefully crafted random noise that evolves over time as a frequency-dependent 
phase shift [Tessendorf02]. Those same algorithms have found their way into 
game code [JensenOl], [Arete03] without significant modification, producing 
beautiful ocean surfaces that evolve at 30 frames per second or more on unexcep¬ 
tional hardware. 

What the FFT algorithms do not give, however, is interactivity between objects 
and the water surface. It would be difficult, for example, to have characters wade 
through a stream and generate a disturbance that depends directly on the motion 
that the player controls. A jet ski thrashing about in the water would not generate 
turbulent waves. Waves in a bathtub cannot bounce back and forth using FFT- 
based simulation. In general, it is not possible in the FFT approach to place an arbi¬ 
trary object in the water and have it interact in a realistic way with the surface 
without a substantial loss of frame rate. For practical purposes, wave surfaces are 
restricted in the ways that the height data can be modified within a frame and 
between frames. 

This article provides a new method, which has been dubbed iWave, for comput¬ 
ing water surface wave propagation that overcomes this limitation. The three scenar¬ 
ios for the stream, jet ski, and bathtub are handled well with iWave. Objects with any 


265 


TeamLRN 



266 


Section 3 Physics 


shape can be present on the water surface and generate waves. Waves that approach 
an object reflect off of it realistically. The entire iWave algorithm amounts to a two- 
dimensional convolution and some masking operations, both suitable for hardware 
acceleration. Even without hardware assistance, a software-only implementation is 
capable of simulating a 128 x 128 water surface height grid at over 30 fps on 1 GHz 
processors. Larger grids will of course slow the frame rate down, and smaller grids 
will speed it up; the speed is directly proportional to the number of grid points. 
Because the method avoids FFTs, it is highly interactive and suitable for a wide range 
of possible applications. 


Linear Waves 

Let’s begin with a quick review of the equations of motion for water surface waves. An 
excellent resource for details on the fluid dynamics is [Kinsman84]. The equations 
that are appropriate here are called the linearized Bernoulli’s equations. The form of 
this equation we use here has a very strange operator that will be explained in a 
moment. The equation is [Tessendorf02] 


d 2 h(x,y, t) Bh(x,y,t) = _ 

"i -v A ' 


dt 


dt 


1 h{x,y,t) 


(3.6.1) 


In this equation, h(x,y, t) is the height of the water surface with respect to the mean 
height at the horizontal position (x, y) at time t. The first term on the left is the vertical 
acceleration of the wave. The second term on the left side, with the constant a, is a 
velocity damping term, not normally a part of the surface wave equation, but is useful 
sometimes to help suppress numerical instabilities that can arise. The term on the right 
side comes from a combination of mass conservation and the gravitational restoring 
force. The operator 



d 2 d 2 

y dx 2 dy 2 


is a mass conservation operator, and we will refer to it as a vertical derivative of the 
surface. Its effect is to conserve the total water mass being displaced. When the height 
of the surface rises in one location, it carries with it a mass of water. To conserve mass, 
there is a region of the surface nearby where the height drops, displacing downward 
the same amount of water that is displaced upward in the first location. 

The next section describes how to evaluate the right-hand side of Equation 3.6.1 
by expressing it as a convolution. Throughout the rest of this article, the height is 
computed on a regular grid, as shown in Figure 3.6.1. The horizontal position {x, y) 
becomes the grid location (i,j) at positions x, = iA and =jA, with the grid spacing A 
the same in both directions. The indices run /=!,..., N and j = 1,. . . , M. 


TeamLRN 



Interactive Water Surfaces 


267 



FIGURE 3.6.1 Layout of the gridfor computing 
wave height. 

Vertical Derivative Operator _ ___ __ 

Like any linear operator acting on a function, the vertical derivative can be imple¬ 
mented as a convolution on the function to which it is applied. In this section, we 
build up this convolution, applied to height data on a regular grid. We also determine 
the best size of the convolution and compute the tap weights. 

As a convolution, the vertical derivative operates on the height grid as 

V-V^ h{i,j) = X X l)h{i + k ’j + l) 

k=-r i=-p (3.6.2) 

The convolution kernel is square, with dimensions (2 P+ 1) X (IP + 1), and can 
be precomputed and stored in a lookup table prior to start of the simulation. The 
choice of the kernel size P affects both the speed and the visual quality of the simula¬ 
tion. The choice P= 6 is the smallest value that gives clearly waterlike motion. 

Figure 3.6.2 shows the kernel elements G(k,0) as a function of k. The two dashed 
vertical lines are at the points k = 6 and k = —6. You can see from the plot that at larger 
values of k , the kernel is mostly zero, and including values of k outside the dashed 
lines would not contribute much to the convolution. If we stop the convolution at 
smaller values, say k = 5 and k = -5, evaluating the convolution is faster, but we will 
miss some small contribution from the k = 6, k = — 6 terms. Experience shows that you 
can get really good-looking waves keeping the terms out to 6, but if you are pressed 


TeamLRN 




268 


Section 3 Physics 


for computation time, stopping the convolution short of that can work, but it will not 
be as visually realistic. Terminating the kernel at a value I K | < 6 sacrifices significant 
amounts of oscillation. This analysis is why the choice P=6 is recommended as the 
best compromise for reasonable wavelike simulation. 


Kernel as a Function of Separation 



FIGURE 3.6.2 The vertical derivative kernel in cross section. Between the two dashed 
lines is the \ K \ < 6 region. 


Computing the kernel values and storing them in a lookup table is a relatively 
straightforward process. The first step is to compute a single number that will scale 
the kernel so that the center value is one. The number is 

G o = ex p( - °#») 

n 

For this sum, q n = ntsq with A q = 0.001 being a good choice for accuracy, and n 
= 1, . . . ,10000. The factor <5 makes the sum converge to a reasonable number, and 
the choice <7=1 works well. With this number in hand, the kernel values are 

G (k,l) = ^ exp(— G^\h(qj) / 


with the parameter r — yk~ + T . The computation time for the kernel elements 
is relatively small, and all of the cost is an initialization; once the elements are com¬ 
puted, they are fixed during the simulation. 


TeamLRN 




Interactive Water Surfaces 


269 


One remaining item needed to compute the convolution kernel elements is a 
formula for the Bessel function J 0 (x). This is included in the C standard math library 
as j 0 ( ). If you do not have access to this, a very convenient approximate fit for this 
function is provided in [Abramowitz72]. Although the formula there is a fitted para¬ 
metric form, it is accurate to within single precision needs, and works well for the 
purposes of this simulation. 

When you perform the convolution at each time step, there are opportunities to 
optimize its speed, both for particular hardware configurations and in software. Soft¬ 
ware optimizations follow because of two symmetries in the convolution kernel: The 
kernel is rotation symmetric, G{k,t) = G(l,k), and the kernel is reflection symmetric 
about both axes, G(k,l) = G(—k,—l) = G(k,—l) = G(-k,l). Without applying any 
symmetries, evaluating the convolution in Equation 3.6.2 directly requires (2 P+ l) 2 
multiplications and additions. Applying these symmetries, the convolution can be 
rewritten as (using the fact that (7(0,0) = 1 by construction) 

h(i + k, j + l) + h{i - k, j - l) + h(i + k,j — l) 

+ h(i-k,j + l) y 

In this form, there are still (2 P+ l) 2 additions, but only (P+ l)P1 2 multiplications. 

This type of convolution can be cast in a form suitable for a SIMD pipeline, so 
graphics cards and DSPs can execute this convolution efficiently. 


j) + X X G ^ k ’ 


k=0 l=t +1 


Wave Propagation _ 

Now that we are able to evaluate the vertical derivative on the height grid, the prop¬ 
agation of the surface can be computed over time. It is simplest to use an explicit 
scheme for time stepping. Although implicit methods can be more accurate and 
stable, they are also slower. Since we are solving a linear equation in this article, an 
explicit approach is fast and stable in the presence of friction, and time step sizes can 
be set to whatever is needed for the display frame rate. If necessary, the friction can be 
kept very low, although for game purposes it might be preferable to have the waves 
dissipate when they are no longer driven by sources. 

To construct the explicit solution, the time derivatives in Equation 3.6.1 must be 
written as finite differences. The second derivative term can be built as a symmetric 
difference, and the dissipative friction term as a forward difference. Rearranging the 
results terms, and assuming a time step A t, the height grid at the next time step is 



Ar) 



2 - aAt 
1 + aAt 



g Af2 y 

l + aAt k ±! p 


p 

£ G(k,l)h(i + 

/=-/> 


1 

1 + aAt 
k,j + /, t) 


(3.6.3) 


TeamLRN 



270 


Section 3 Physics 


In terms of data structures, this algorithm for propagation can be run with three 
copies of the heightfield grid. For this discussion, the grids are held in the float arrays 
height, vertical_derivative, and previous_height. During the simulation, the 
array height always holds the new height grid, previous_height holds the height 
grid from the previous time step, and vertical_derivative holds the vertical deriv¬ 
ative of the height grid from the previous time step. Before simulation begins, they 
should all been initialized to zero for each element. The pseudocode to accomplish 
the propagation is 

float height[N*M]; 

float vertical_derivative[N*M]; 

float previous_height[N*M]; 

// ... initialize to zero ... 

// ... begin loop over frames ... 

// --- This is the propagation code --- 
// Convolve height with the kernel 
// and put it into vertical_derivative 
Convolve( height, vertical_derivative ); 

float temp; 

for(int k=0; k<N*M; k++) 

{ 

temp = height[k]; 

height[k] = height[k]*(2.0- 

alpha*dt)/(1.0+alpha*dt) 

- previous_height[k]/(1.0+alpha*dt) 

- vertical_derivative[k] 

*g*dt*dt/(1.0+alpha*dt); 

previous_height[k] = temp; 

} 

// — end propagation code --- 
// ... end loop over frames ... 

The quantities in vertical_derivative and previous_height could be useful for 
embellishing the visual look of the waves. For example, a large value in 
vertical_derivative indicates strong gravitational attraction of the waves back to 
the mean position. Comparing the value in previous_height with that in height at 
the location of strong vertical_derivative can determine roughly whether the wave 
is at a peak or a trough. If it is at a peak, a foam texture could be used in that area. 
This is not a concrete algorithm grounded in physics or oceanography, but just a spec¬ 
ulation about how peaks of the waves might be found. The point of this is simply that 
the two additional grids vertical_derivative and previous_height could have some 
additional benefit in the simulation and rendering of the wave height field beyond 
just the propagation steps. 


TeamLRN 




Interactive Water Surfaces 


271 


Interacting Obstructions and Sources 

Up to this point, we have built a method to propagate waves in a water surface simu¬ 
lation. While the propagation involves a relatively fast convolution, everything we 
have discussed could have been accomplished just as efficiently (possibly more effi¬ 
ciently) with a FFT approach such as the ones mentioned in the introduction. The 
real power of this convolution method is the ease with which some additional 2D 
processing can generate highly realistic interactions between objects in the water, and 
pump disturbances into the water surface. 

The fact that we can get away with 2D processing to produce interactivity is, in 
some ways, a miracle. Normally in a fluid dynamic simulation the fluid velocity on and 
near a boundary is reset according to the type of boundary condition and requires 
understanding of geometric information about the boundary, such as its outward nor¬ 
mal. Flere we get away with effectively none of that analysis, which is critical to the 
speed of this approach. 

Sources 

One way of creating motion in the fluid is to have sources of displacement. A source 
is represented as a 2D grid s(i,j) of the same size and dimensions as the height grid. 
The source grid should have zero values wherever no additional motion is desired. At 
locations in which the waves are being “poked” and/or “pulled,” the value of the 
source grid can be positive or negative. Then, just prior to propagation step in Equa¬ 
tion 3.6.3, the height grid is updated with h(i,j) = h(i,j) + Since the source is an 
energy input per frame, it should change over the course of the simulation, unless a 
constant buildup of energy really is desired. An impulse source generates a ripple. 

Obstructions 

Obstructions are extremely easy to implement in this scheme. An additional grid for 
obstructions is filled with float values, primarily with two extreme values. This grid 
acts as a mask indicating where obstructions are present. At each grid point, if there is 
no obstruction present, then the value of the obstruction grid at that point is 1.0. If a 
grid point is occupied by an obstruction, then the obstruction grid value is 0.0. At 
grid points on the border around an obstruction, the value of the obstruction grid is 
some intermediate value between 0.0 and 1.0. The intermediate region acts as an 
anti-aliasing of the edge of the obstruction. 

Given this obstruction mask, the obstructions influence is computed by simply 
multiplying the height grid by the obstruction mask, so that the wave height is 
forced to zero in the presence of the obstruction, and left unchanged in areas outside 
the obstruction. Amazingly, that is all that must be done to properly account for 
objects on the water surface! This simple step causes waves that propagate to the 
obstruction to reflect correctly off it. It also produces refraction of waves that pass 


TeamLRN 



272 


Section 3 Physics 


through a narrow slit channel in an obstruction. In addition, it permits the obstruc¬ 
tion to have any shape at all, animating in any way that the user wants it to. 

The pseudocode for an application with sources and obstructions is: 

float source[N*M], obstruction[N*M]; 

// ... set the source and obstruction grids 

for(int k = 0; k < N*M; k++) 

{ 

height[k] += source[k]; 

height[k] *= obstruction^]; 

} 

II ... now apply propagation 

Wakes 

Wakes from moving objects are naturally produced by the iWave method of interactiv¬ 
ity. In this special case, the shape of the obstruction is also the shape of the source. Set¬ 
ting source[k] = 1.0-obstruction[k] works as long as there is an anti-aliased region 
around the edge of the obstruction. With this choice, moving an obstacle around 
in the grid produces a wake behind it that includes the V-shaped Kelvin wake. It also 
produces a type of stern wave and waves running along the side of the obstacle. The 
details of the shape, timing, and extent of these wake components are sensitive to 
the shape and motion of the obstacle. 

Ambient Waves 

The iWave method is not very effective at generating persistent large-scale wave phe¬ 
nomena like open ocean waves. If an application desires these “ambient waves” that 
are not generated in the iWave method, there is an additional procedure that avoids 
explicitly simulating the ambient waves. 

The ambient waves consist of a height grid that has been generated by some other 
procedure. For example, FFT methods could be used to generate ocean waves and store 
them in a height grid. Since we are only trying to compute the interaction of the ambi¬ 
ent waves with an obstruction, the ambient waves should not contribute to the simula¬ 
tion outside the region of the obstruction. The pseudocode for modifying the height 
grid, prior to propagation and just after application of obstructions and sources as done 
previously, is 

float ambient[N*M]; 

// ... set the ambient grid for this time step 

// ... just after the source and obstruction, apply: 

for (int k = 0; k < N*M; k++) 

{ 

height[k] -= ambient[k]*(1.0-obstruction[k]); 

} 


TeamLRN 



Interactive Water Surfaces 


273 


// ... now apply the propagation 

With this method, ambient waves of any character can interact with objects of 
any animating shape. 

Grid Boundaries 

Up to this point, we have ignored the problem of how to treat the boundaries of 
the grid. The problem is that the convolution kernel requires data from grid points 
a distance P in all four directions from the central grid point of the convolution. 
Therefore, when the central grid point is fewer than P points from a boundary of 
the grid, missing data must be fdled in according to some criterion. There 
are two types of boundary conditions that are fairly easy to apply: periodic and 
reflecting boundaries. 

Periodic Boundaries 

In this situation, a wave encountering a boundary appears to continue to propagate 
inward from the boundary on the opposite side. In performing the convolution near 
the boundaries, the grid coordinates in Equation 3.6.3 (z + k and j + l) may be outside 
of the ranges [0,7V- 1] and [0 ,M - 1]. Applying the modulus (i+k)%N is guaranteed 
to be in the range [—TV 4- 1,7V— 1]. To insure that the result is always positive, a dou¬ 
ble modulus can be used: ( (i+k)%N + N)%N. Doing the same for they + /coordinate 
insures that periodic boundary conditions are enforced. 

Reflecting Boundaries 

Reflecting boundaries turn a wave around and send it back into the grid from the bound¬ 
ary the wave is incident on, much like a wave that reflects off an obstacle in the water. If 
the coordinate i + k is greater than TV— 1 then it is changed to 27V—z—A If the coordinate 
is less than zero, it is negated; in other words, it becomes —i—k, which is positive. An iden¬ 
tical procedure should also be applied to they + / coordinate. 

To efficiently implement either of these two types of boundary treatments, the 
fastest approach is to divide the grid into nine regions as follows: 

1. The inner portion of the grid with the range of coordinates ze [P,N-\-P] 
and ye [P,M-\-P\. 

2. The right-hand side ze [7V-.P,7V-1] and ye [P,M—\—P\. 

3. The left-hand side ze [0,/*—1] and ye [P,M—\—P], 

4. The top side ze [P,N-\-P] and ye [0,P-1]. 

5. The bottom side ze [P,7V— 1— P\ and ye [M—P,M— 1 ]. 

6. The four corners that remain. 

Within each region, the particular boundary treatment required can be coded 
efficiently without conditionals or extra modulus operations. 


TeamLRN 



274 


Section 3 Physics 


Surface Tension 

So far, the type of simulation we discussed is the propagation of gravity waves. Grav¬ 
ity waves dominate surface flows on scales of approximately a foot or larger. On 
smaller scales, the character of the propagation changes to include surface tension. 
Surface tension causes waves to propagate faster at smaller spatial scales, which tends 
to make the surface appear to be more rigid than without it. For our purposes, surface 
tension is characterized by a length scale L T , which determines the maximum size of 
the surface tension waves. The only change required of our procedure is a different 
computation of the convolution kernel. The kernel calculation becomes 

G ( k > l) = X £ V 1 + <7*4 ex P )jo{ c Jn r ) / G o 
n / 

Other than this change, the entire iWave process is the same. 


Conclusion 


The iWave method of water surface propagation is a very flexible approach to creating 
interactive disturbances of water surfaces. Because it is based on 2D convolution and 
some simple 2D image manipulation, high frame rates can be obtained even in a 
software-only implementation. Hardware acceleration of the convolution should 
make iWave suitable for many game platforms. The increased interactivity of the 
water surface with objects in a game could open new areas of gameplay that previously 
were not available to the game developer. 

References 

[Abramowitz72] Abramowitz, Milton, and Irene A. Stegun, Handbook of Mathemati¬ 
cal Functions, Dover, 1972. Sections 9.4.1 and 9.4.3. 

[Arete03] Arete Entertainment, available online at www.areteis.com. 

[JensenOl] Jensen, Lasse, online tutorial, available online at www.gamasutra.com/gdce/ 
jensen/jensen_01 -htm, 2001. 

[Kinsman84] Kinsman, Blair, Wind Waves, Dover, 1984. 

[Tessendorf02] Tessendorf, Jerry, “Simulating Ocean Water,” Simulating Nature, SIG- 
GRAPH Course Notes, available online at http://users, adelphia. net !- tessendorf, 2002. 


TeamLRN 



*/