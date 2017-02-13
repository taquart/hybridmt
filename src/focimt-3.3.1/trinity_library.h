/*
 * trinity_library.h
 *
 *  Created on: 8 cze 2016
 *      Author: gregus
 */

#ifndef TRINITY_LIBRARY_H_
#define TRINITY_LIBRARY_H_

#include <string>
#include <algorithm>
#include <stdio.h>
#include <math.h>

//=============================================================================
//=============================================================================
//=============================================================================

namespace Taquart {
  /*! \defgroup stdexception STD-compliant exception classes (libexceptions.a)
   *  \ingroup trilib
   */

//! Parent class for simple exception classes.
  /*! TriException class is a parent class for various exception classes
   *  in Trinity namespace.
   *
   *  \version 1.4.0 [2006.10.01] Removed <stdio.h> header file, TriException
   *   class was totally rebuilt (copy constructor and assignment operator
   *   added).
   *  \version 1.2.0 [2006.09.16] Removed <string> header file.
   *  \version 1.1.2 [2006.06.13] Corrected header file.
   *  \version 1.1.1 [2005.11.06] Fixed a few errors in the code, full
   *   documentation included.
   *  \version 1.1.0 [2005.11.05] Re-coded from scratch
   *  \version 1.0.0 [2004.11.10] First version released
   *
   *  \ingroup stdexception
   *  \ingroup trilib
   */
  class TriException {
    public:
      //! Pointer to the exception message.
      char * Message;

      //! Default constructor.
      TriException(void);

      //! Constructor.
      /*! \param AMessage Custom exception message (will be merged with
       *  a default error message for parent class.)
       */
      TriException(const char * AMessage);

      //! Copy constructor.
      /*! \param ASource Reference to the source object.
       */
      TriException(const TriException& ASource);

      //! Assignment operator.
      /*! \param ASource Reference to the source object.
       *  \return Reference to the current object.
       */
      TriException& operator=(const TriException& ASource);

      //! Default destructor.
      virtual ~TriException(void);

    private:

      unsigned int length(const char * Text);
      void clear(void);
      void fill(const char * AMessage);

    protected:
  };

//! "Wrong operation" exception class.
  /*! Trinity::TriEWrongOperation is thrown when forbidden operation occurs.
   *  \ingroup stdexception
   *  \ingroup trilib
   */
  class TriEWrongOperation: public TriException {
    public:
      //! Default constructor
      /*! \param AMessage Exception message.
       */
      TriEWrongOperation(const char * AMessage = "Wrong operation");

    private:
    protected:
  };

//! "Out of range" exception class.
  /*! Trinity::TriEOutOfRange is thrown when the number is out of range.
   *  \ingroup stdexception
   *  \ingroup trilib
   */
  class TriEOutOfRange: public Taquart::TriException {
    public:
      //! Default constructor
      /*! \param AMessage Exception message.
       */
      TriEOutOfRange(const char * AMessage = "Out of range");

    private:
    protected:
  };

//! "Operation aborted" exception class.
  /*! Trinity::TriEOperationAborted is thrown when current operation has
   *  been aborted.
   *  \ingroup stdexception
   *  \ingroup trilib
   */
  class TriEOperationAborted: public TriException {
    public:
      //! Default constructor
      /*! \param AMessage Exception message.
       */
      TriEOperationAborted(const char * AMessage = "Operation aborted");

    private:
    protected:
  };

//! "Input/output" exception class.
  /*! Trinity::TriEIOError is thrown when a general input/output error occurs.
   *  \ingroup stdexception
   *  \ingroup trilib
   */
  class TriEIOError: public TriException {
    public:
      //! Default constructor
      /*! \param AMessage Exception message.
       */
      TriEIOError(const char * AMessage = "Input / Output error");

    private:
    protected:
  };

//! "Conversion failed" exception class.
  /*! Trinity::TriEConversionError is thrown when conversion from a one type to
   *  another one is not available (e.g. while converting string to
   *  floating-point number is impossible).
   *  \ingroup stdexception
   *  \ingroup trilib
   */
  class TriEConversionError: public TriException {
    public:
      //! Default constructor
      /*! \param AMessage Exception message.
       */
      TriEConversionError(const char * AMessage = "Conversion error");

    private:
    protected:
  };

//! "Empty string" exception class.
  /*! Trinity::TriEEmptyString is thrown when some string is empty (e.g. its
   *  length is equal to zero). To signal that the pointer is NULL, use
   *  Trinity::TriENullPointer instead.
   *  \ingroup stdexception
   *  \ingroup trilib
   */
  class TriEEmptyString: public TriException {
    public:
      //! Default constructor
      /*! \param AMessage Exception message.
       */
      TriEEmptyString(const char * AMessage = "Empty string");

    private:
    protected:
  };

//! "NULL pointer" exception class.
  /*! Trinity::TriENullPointer is thrown when parameter specified in e.g.
   *  function call is NULL. To signal that the string is empty, use
   *  Trinity::TriENullPointer instead.
   *  \ingroup stdexception
   *  \ingroup trilib
   */
  class TriENullPointer: public TriException {
    public:
      //! Default constructor
      /*! \param AMessage Exception message.
       */
      TriENullPointer(const char * AMessage = "NULL Pointer");

    private:
    protected:
  };
} // namespace Trinity


//=============================================================================
//=============================================================================
//=============================================================================

namespace Taquart {

  class String {
    public:
      String(void);
      String(const std::string &AText);
      String(const char *AText);
      String(const Taquart::String &Source);
      String(const double Source);
      Taquart::String &operator= (const Taquart::String &Source);
      Taquart::String operator+ (const Taquart::String &Source);
      char &operator[] (int i);
      void Assign(const Taquart::String &Source);
      virtual ~String(void);
      const char * c_str(void);
      std::string c_txt(void);
      int Pos(const std::string &Needle);
      int Pos(const char *Needle);
      int Pos(Taquart::String &Needle);
      //int Pos(char Needle);
      int Length(void);
      Taquart::String SubString(int Pos, int Length);
      Taquart::String Trim(void);
      Taquart::String TrimLeft(void);
      Taquart::String TrimRight(void);
      Taquart::String UpperCase(void);
      Taquart::String LowerCase(void);
      double ToDouble(void);
      int ToInt(void);
      //friend bool operator!= (Taquart::String &Object1, Taquart::String &Object2);
      friend bool operator== (const Taquart::String &Object1, const Taquart::String &Object2);
    private:
      std::string Text;
    protected:
  };

  Taquart::String ExtractFileName(Taquart::String Filename);

  Taquart::String FormatFloat(Taquart::String Format, double Value);
  //Taquart::String FormatFloat(Taquart::String Format, int Value);

  bool operator== (const Taquart::String &Object1, const Taquart::String &Object2);
} /* namespace Taquart */

//=============================================================================
//=============================================================================
//=============================================================================

namespace Taquart {
  double alog(double x);
  double alog10(double x);
  double amax1(double x, double y);
  double amin1(double x, double y);
  double sign(double x, double y);
  double datan2(double y, double x);
}

//=============================================================================
//=============================================================================
//=============================================================================

namespace Taquart {
  //! Returns the sample standard deviation for elements in an array.
  /*! Trinity::std calculates the sample standard deviation (the square root
   *  of the sample variance) of all values in the \p X array parameter.
   *  \p Size indicates the number of elements in the array.
   *  \param X Pointer to the array.
   *  \param Size Array size.
   *  \return Standard deviation value.
   *  \ingroup trilib
   */
  double std(double * X, unsigned int Size) throw (Taquart::TriENullPointer,
      Taquart::TriEOutOfRange);

  //! Returns the average of all values in an array.
  /*! Trinity::mean calculates the arithmetic average of all the values in
   *  the \p X array parameter. The \p Size parameter gives the number
   *  of array elements.
   *  \param X Pointer to the array.
   *  \param Size Array size.
   *  \return Mean value.
   *  \ingroup trilib
   */
  double mean(double * X, unsigned int Size) throw (Taquart::TriENullPointer,
      Taquart::TriEOutOfRange);
}

//=============================================================================
//=============================================================================
//=============================================================================

namespace Taquart {
  //! Color description class.
  /*! \ingroup tricairo
   */
  typedef class TriCairo_Color {
    public:
      double R;
      double G;
      double B;
      double A;
      TriCairo_Color(void);
      TriCairo_Color(double R, double G, double B, double A = 1.0);
      TriCairo_Color(unsigned char R, unsigned char G, unsigned char B,
          unsigned char A = 255);
      //TriCairo_Color(TColor t, double A=1.0);
      void Dispatch(double &r, double &g, double &b, double &a);
    private:
    protected:
  } TCColor;
}

//=============================================================================
//=============================================================================
//=============================================================================

#include <cairo/cairo.h>
//#include <cairo/cairo-win32.h>
#include <cairo/cairo-svg.h>
#include <cairo/cairo-pdf.h>
#include <cairo/cairo-ps.h>

//! \defgroup tricairo Interface to CAIRO library.

namespace Taquart {

  //! Type of line joint
  /*! \ingroup tricairo
   */
  enum TriCairo_LineJoin {
    ljMiter, ljBevel, ljRound
  };

  //! Type of line cap.
  /*! \ingroup tricairo
   */
  enum TriCairo_LineCap {
    lcButt, lcRound, lcSquare
  };

  //! Type of output canvas.
  /*! \ingroup tricairo
   */
  enum TriCairo_CanvasType {
    ctBitmap, ctSurface, ctSVG, ctPDF, ctPS
  };

  //! Font style.
  /*! \ingroup tricairo
   */
  enum TriCairo_FontStyle {
    fsNormal, fsItalic, fsOblique, fsNormalBold, fsItalicBold, fsObliqueBold
  };

  //! Horizontal alignment of text.
  /*! \ingroup tricairo
   */
  enum TriCairo_HorizontalAlignment {
    haLeft, haCenter, haRight
  };

  //! Vertical alignment of text.
  /*! \ingroup tricairo
   */
  enum TriCairo_VerticalAlignment {
    vaTop, vaMiddle, vaBottom
  };

  //! Base class wrapping the interface between BDS2006 and Cairo library.
  /*! \ingroup tricairo
   */
  class TriCairo {
    public:
      TriCairo(unsigned int width, unsigned int height,
          Taquart::TriCairo_CanvasType canvastype,
          Taquart::String Filename = "");
      virtual ~TriCairo(void);
      virtual void Save(Taquart::String filename);
      //Graphics::TBitmap * GetBitmap(void);
      //Graphics::TBitmap * CreateBitmap(void);
      //void DrawCanvas(TCanvas * Canvas, int Left = 0, int Top = 0);

      // Drawing functions.
      void Color(Taquart::TriCairo_Color Color);
      void ColorD(double r, double g, double b, double a = 1.0);
      void ColorB(unsigned int r, unsigned int g, unsigned int b,
          unsigned int a = 255);
      void Clear(double r, double g, double b, double a = 1.0);
      void Clear(Taquart::TriCairo_Color Color);
      void MoveTo(double x, double y);
      void LineTo(double x, double y);
      void LineToRel(double x, double y);
      void LineWidth(double w);
      void LineCap(TriCairo_LineCap lc);
      void LineJoin(TriCairo_LineJoin lj);
      void Font(Taquart::String Name, double Size, TriCairo_FontStyle fs);
      void Text(double x, double y, Taquart::String Text,
          TriCairo_HorizontalAlignment ha = haLeft,
          TriCairo_VerticalAlignment va = vaTop);
      void Stroke(void);
      void StrokePreserve(void);
      void Fill(void);
      void FillPreserve(void);
      void Arc(double x, double y, double r, double start = 0.0,
          double end = 2 * M_PI);
      void Rectangle(double x, double y, double w, double h);
      void Circle(double x, double y, double r);
      void ClosePath(void);

      Taquart::String Filename;
    private:
      //Graphics::TBitmap * Bitmap;

    protected:
      cairo_surface_t * CreateSurface(TriCairo_CanvasType CanvasType);

      // Additional functions.
      TriCairo_CanvasType CanvasType;
      const unsigned int Width;
      const unsigned int Height;
      cairo_surface_t * surface;
      cairo_t * cr;
  };

}

//=============================================================================
//=============================================================================
//=============================================================================

#define EPSIL 0.0001
#define DEG2RAD (M_PI / 180.0)
#define sind(x) sin ((x) * DEG2RAD)
#define cosd(x) cos ((x) * DEG2RAD)
#define tand(x) tan ((x) * DEG2RAD)

namespace Taquart {

  typedef struct nodal_plane {
      double str;
      double dip;
      double rake;
  } TriCairo_NodalPlane;

  typedef struct TriCairo_Axis {
      double str;
      double dip;
      double val;
      TriCairo_Axis(void) {
        str = 0.0;
        dip = 0.0;
        val = 0.0;
      }
  } AXIS;

  struct MOMENT {
      double mant;
      int exponent;
  };

  struct MECHANISM {
      struct nodal_plane NP1;
      struct nodal_plane NP2;
      struct MOMENT moment;
      double magms;
  };

  typedef struct MECHANISM st_me;

  void axe2dc(AXIS T, AXIS P, nodal_plane *NP1,
   nodal_plane *NP2);
  void StrikeDipRake2MT(double S, double D, double R, double &M11, double &M22,
      double &M33, double &M12, double &M13, double &M23);
  double zero_360(double str);
  void sincosd(double i, double *s, double *c);
  void sincos(double i, double *s, double *c);
  double computed_dip1(nodal_plane NP1);
  double computed_rake1(struct nodal_plane NP1);
  double computed_strike1(struct nodal_plane NP1);
  void define_second_plane(struct nodal_plane NP1, struct nodal_plane *NP2);
  double null_axis_dip(double str1, double dip1, double str2, double dip2);
  double null_axis_strike(double str1, double dip1, double str2, double dip2);
  double computed_rake2(double str1, double dip1, double str2, double dip2,
      double fault);
}
//=============================================================================
//=============================================================================
//=============================================================================

namespace Taquart {
  typedef int GMT_LONG;

  //! Hemisphere projection.
  /*! \ingroup tricairo
   */
  enum TriCairo_Hemisphere {
    heLower = 0, heUpper = 1
  };

  //! Type of station marker.
  /*! \ingroup tricairo
   */
  enum TriCairo_MarkerType {
    mtCircle = 0, mtSquare = 1, mtPlusMinus = 2, mtBWCircle = 3
  };

  //! Type of network projection.
  /*! \ingroup tricairo
   */
  enum TriCairo_Projection {
    prWulff = 0, prSchmidt = 1
  };

  //! Structure stores information about plunge and trend of axis.
  /*! \ingroup tricairo
   */

  //! Structure stores information about strike, dip and rake of a fault.
  /*! \ingroup tricairo
   */

  /*
   typedef struct DLL_EXP TriCairo_NodalPlane
   {
   double str;
   double dip;
   double rake;
   } nodal_plane;
   */
  //! Structure stores information moment tensor.
  /*! The corresponding matrix elements correspond to the moment tensor
   *  components according to CMT convention.
   *  \ingroup tricairo
   */
  typedef struct TriCairo_MomentTensor {
      double f[6]; /* mrr mtt mff mrt mrf mtf in 10**expo dynes-cm */
      TriCairo_MomentTensor(void) {
        for (int i = 0; i < 6; i++)
          f[i] = 0;
      }
  } M_TENSOR;

  //! Class for producing the graphical representation of moment tensor component.
  /*! This class is capable to produce a graphical representation of the seismic
   *   moment tensor (so called beach balls).
   *  \ingroup tricairo
   */
  class TriCairo_Meca: public TriCairo {
    public:
      // Constructor.
      TriCairo_Meca(unsigned int width, unsigned int height,
          TriCairo_CanvasType type, Taquart::String filename = "");
      // Destructor.
      virtual ~TriCairo_Meca(void);
      //void Draw(double M11, double M12, double M13, double M22, double M23, double M33,
      //  double s1, double d1, double r1, double s2, double d2, double r2);

      // Public variables.
      bool DrawCross;
      bool DrawDC;

      TriCairo_Hemisphere Hemisphere;
      unsigned int Margin; // Margin size.
      unsigned int BWidth;
      unsigned int BHeight;
      unsigned int BRadius;
      unsigned int BXo;
      unsigned int BYo;
      double BOutlineWidth;
      TCColor BOutlineColor;
      TCColor BPlusColor; // For compressional part.
      TCColor BMinusColor; // For dilatational part.
      TCColor BTensorOutline;
      TCColor BDCColor;
      double BDCWidth;

      bool DrawAxis;
      double AxisFontSize;
      Taquart::String AxisFontFace;

      bool DrawStations;
      bool DrawStationName;
      double StationMarkerSize;
      double StationFontSize;
      TCColor StationPlusColor;
      TCColor StationMinusColor;
      TriCairo_MarkerType StationMarkerType;
      Taquart::String StationFontFace;
      TCColor StationTextColor;

      TriCairo_Projection Projection;

      // Main drawing routines.
      void Station(double Azimuth, double Takeoff, double Disp,
          Taquart::String Label, double &mx, double &my, double error = 0.0);
      void Tensor(AXIS T, AXIS N, AXIS P);
      void Axis(AXIS A, Taquart::String Text);
      void DoubleCouple(double Strike, double Dip);
      void DrawStationMarker(double x, double y, double disp,
          Taquart::String Label);
      void CenterCross(void);
      void GMT_momten2axe(M_TENSOR mt, AXIS *T, AXIS *N, AXIS *P);
      void Station(double GA[], double Disp, Taquart::String Label, double &mx,
          double &my, double error = 0.0);

    protected:

      // Upper or lower hemisphere projection routines.
      void Project(double &X, double &Y);
      void Project(double * X, double * Y, unsigned int npoints);

    private:
      // Additional routines.
      double squared(double v);
      void axe2dc(AXIS T, AXIS P, nodal_plane *NP1, nodal_plane *NP2);
      double proj_radius2(double str1, double dip1, double str);
      void ps_circle(double x0, double y0, double radius_size, TCColor fc);
      void Polygon(double xp1[], double yp1[], int npoints, TCColor oc,
          bool fill, TCColor fc = TCColor(), double OutlineWidth = 1.0);

      GMT_LONG GMT_jacobi(double *a, GMT_LONG *n, GMT_LONG *m, double *d,
          double *v, double *b, double *z, GMT_LONG *nrots);

      void * GMT_memory(void *prev_addr, GMT_LONG nelem, size_t size);
      void GMT_free(void *addr);

  };
// class TriCairo_Meca
}

//=============================================================================
//=============================================================================
//=============================================================================

//=============================================================================
//=============================================================================
//=============================================================================

#endif /* TRINITY_LIBRARY_H_ */
