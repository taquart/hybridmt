/*
 * trinity_library.cpp
 *
 *  Created on: 8 cze 2016
 *      Author: gregus
 */

#include <stdlib.h>
#include "trinity_library.h"
using namespace Taquart;

//=============================================================================
//=============================================================================
//=============================================================================

TriException::TriException(void) {
  Message = 0;
}

//---------------------------------------------------------------------------
TriException::TriException(const char * AMessage) {
  Message = 0;
  fill(AMessage);
}

//---------------------------------------------------------------------------
TriException::TriException(const TriException& ASource) {
  Message = 0;
  fill(ASource.Message);
}

//---------------------------------------------------------------------------
TriException& TriException::operator=(const TriException& ASource) {
  if (this != &ASource) {
    fill(ASource.Message);
  }
  return *this;
}

//---------------------------------------------------------------------------
TriException::~TriException(void) {
  clear();
}

//---------------------------------------------------------------------------
unsigned int TriException::length(const char * Text) {
  unsigned int i = 0;
  while (*(Text + i) != 0)
    i++;
  return i;
}

//---------------------------------------------------------------------------
void TriException::clear(void) {
  if (Message) {
    delete Message;
    Message = 0;
  }
}

//---------------------------------------------------------------------------
void TriException::fill(const char * AMessage) {
  clear();
  int Len = length(AMessage);
  Message = new char[Len + 1];
  for (int i = 0; i < Len; i++) {
    Message[i] = AMessage[i];
  }
  Message[Len] = 0;
}

//---------------------------------------------------------------------------
TriEWrongOperation::TriEWrongOperation(const char * AMessage) :
    TriException(AMessage) {
}

//---------------------------------------------------------------------------
TriEOutOfRange::TriEOutOfRange(const char * AMessage) :
    TriException(AMessage) {
}

//---------------------------------------------------------------------------
TriEOperationAborted::TriEOperationAborted(const char * AMessage) :
    TriException(AMessage) {
}

//---------------------------------------------------------------------------
TriEIOError::TriEIOError(const char * AMessage) :
    TriException(AMessage) {
}

//---------------------------------------------------------------------------
TriEConversionError::TriEConversionError(const char * AMessage) :
    TriException(AMessage) {
}

//---------------------------------------------------------------------------
TriEEmptyString::TriEEmptyString(const char * AMessage) :
    TriException(AMessage) {
}

//---------------------------------------------------------------------------
TriENullPointer::TriENullPointer(const char * AMessage) :
    TriException(AMessage) {
}

//=============================================================================
//=============================================================================
//=============================================================================

Taquart::String Taquart::FormatFloat(Taquart::String Format, double Value) {
  char buffer[250];
  snprintf(buffer, 249, Format.c_str(), Value);
  return Taquart::String(buffer);
}

//-----------------------------------------------------------------------------
Taquart::String Taquart::ExtractFileName(Taquart::String Filename) {
  // Remove directory if present.
  // Do this before extension removal incase directory has a period character.
  std::string filename(Filename.c_str());
  const size_t last_slash_idx = filename.find_last_of("\\/");
  if (std::string::npos != last_slash_idx) {
    filename.erase(0, last_slash_idx + 1);
  }

  // Remove extension if present.
  const size_t period_idx = filename.rfind('.');
  if (std::string::npos != period_idx) {
    filename.erase(period_idx);
  }
  return Taquart::String(filename);
}

//-----------------------------------------------------------------------------
Taquart::String::String(void) {
  // Default constructor.
}

Taquart::String Taquart::String::UpperCase(void) {
  std::string Out(Text);
  std::transform(Out.begin(), Out.end(), Out.begin(), ::toupper);
  return Out;
}

Taquart::String Taquart::String::LowerCase(void) {
  std::string Out(Text);
  std::transform(Out.begin(), Out.end(), Out.begin(), ::tolower);
  return Out;
}

Taquart::String::String(const std::string &AText) {
  Text = AText;
}

Taquart::String::String(const char *AText) {
  Text = std::string(AText);
}

std::string Taquart::String::c_txt(void) {
  return Text;
}

Taquart::String::String(const double Source) {
  Text = FormatFloat("%f", Source).c_txt();
}

const char * Taquart::String::c_str(void) {
  return Text.c_str();
}

int Taquart::String::Length(void) {
  return Text.length();
}

int Taquart::String::Pos(const char *Needle) {
  std::string a(Needle);
  return Pos(a);
}

/*
 int Taquart::String::Pos(char character) {
 std::string a(&character);
 return Pos(a);
 }
 */

double Taquart::String::ToDouble(void) {
  return atof(Text.c_str());
}

int Taquart::String::ToInt(void) {
  return atoi(Text.c_str());
}

int Taquart::String::Pos(Taquart::String &Needle) {
  return Pos(Needle.c_str());
}

int Taquart::String::Pos(const std::string &Needle) {
  std::size_t position = Text.find(Needle);
  if (position != std::string::npos)
    return position + 1;
  else
    return 0;
}

Taquart::String::String(const Taquart::String &Source) {
  Assign(Source);
}

Taquart::String& Taquart::String::operator=(const Taquart::String &Source) {
  Assign(Source);
  return *this;
}

Taquart::String Taquart::String::operator+(const Taquart::String &Source) {
  return Taquart::String(Text + Source.Text);
}

void Taquart::String::Assign(const Taquart::String &Source) {
  Text = Source.Text;
}

Taquart::String Taquart::String::SubString(int pos, int len) {
  return Taquart::String(Text.substr(pos - 1, len).c_str());
}

Taquart::String::~String(void) {
  // Default destructor.
}

Taquart::String Taquart::String::Trim(void) {
  return TrimRight().TrimLeft();
}

Taquart::String Taquart::String::TrimLeft(void) {
  size_t startpos = Text.find_first_not_of(" \n\r\t");
  if (startpos == std::string::npos)
    return Taquart::String();
  else
    return Taquart::String(Text.substr(startpos));
}

Taquart::String Taquart::String::TrimRight(void) {
  size_t endpos = Text.find_last_not_of(" \n\r\t");
  return
      (endpos == std::string::npos) ?
          Taquart::String() : Taquart::String(Text.substr(0, endpos + 1));
}

char & Taquart::String::operator[](int i) {
  return Text[i - 1];
}

bool Taquart::operator==(const Taquart::String &Object1,
    const Taquart::String &Object2) {
  return Object1.Text == Object2.Text;
}

//=============================================================================
//=============================================================================
//=============================================================================

double Taquart::alog(double Value) {
  return log(Value);
}

//-----------------------------------------------------------------------------
double Taquart::alog10(double Value) {
  return log10(Value);
}

//-----------------------------------------------------------------------------
double Taquart::amax1(double value1, double value2) {
  return value1 > value2 ? value1 : value2;
}

//-----------------------------------------------------------------------------
double Taquart::amin1(double value1, double value2) {
  return value1 > value2 ? value2 : value1;
}

//-----------------------------------------------------------------------------
double Taquart::sign(double value1, double value2) {
  if (value2 >= 0.0)
    return fabs(value1);
  else
    return -fabs(value1);
}

//-----------------------------------------------------------------------------
double Taquart::datan2(double y, double x) {
  double arctg;
  const double rdeg = 180.0 / M_PI;

  if (fabs(x) < 0.0001) {
    if (fabs(y) < 0.0001) {
      arctg = 0.0;
    }
    else
      arctg = y < 0.0 ? -90.0 : 90.0;
  }
  else if (x < 0.0)
    arctg = y < 0.0 ? atan(y / x) * rdeg - 180.0 : atan(y / x) * rdeg + 180.0;
  else
    arctg = atan(y / x) * rdeg;

  return (arctg);
}

//=============================================================================
//=============================================================================
//=============================================================================

double Taquart::mean(double * X, unsigned int Size)
    throw (Taquart::TriENullPointer, Taquart::TriEOutOfRange) {
  if (X == 0L)
    throw Taquart::TriENullPointer("Trinity:mean(): Input pointer is NULL.");
  else if (Size == 0)
    throw Taquart::TriEOutOfRange(
        "Trinity:mean(): Size of input vector is zero.");
  else {
    double s = 0.0f;
    for (unsigned int i = 0; i < Size; i++)
      s += *(X + i);
    return s / double(Size);
  }
}

//---------------------------------------------------------------------------
double Taquart::std(double * X, unsigned int Size)
    throw (Taquart::TriENullPointer, Taquart::TriEOutOfRange) {
  if (X == 0L)
    throw Taquart::TriENullPointer("Trinity:std(): Input pointer is NULL.");
  else if (Size == 0)
    throw Taquart::TriEOutOfRange(
        "Trinity:std(): Size of input vector is zero.");
  else if (Size == 1)
    return 0.0f;
  else {
    const double MeanValue = mean(X, Size);
    double s = 0.0f;
    for (unsigned int i = 0; i < Size; i++)
      s += ((*(X + i) - MeanValue) * (*(X + i) - MeanValue));
    return sqrt(s / double(Size - 1));
  }
}
//=============================================================================
//=============================================================================
//=============================================================================

TriCairo_Color::TriCairo_Color(void) {
  R = 0.0;
  G = 0.0;
  B = 0.0;
  A = 1.0;
}

//---------------------------------------------------------------------------
TriCairo_Color::TriCairo_Color(double r, double g, double b, double a) {
  R = r;
  G = g;
  B = b;
  A = a;
}

//---------------------------------------------------------------------------
TriCairo_Color::TriCairo_Color(unsigned char r, unsigned char g,
    unsigned char b, unsigned char a) {
  R = (double) r / 255.0;
  G = (double) g / 255.0;
  B = (double) b / 255.0;
  A = (double) a / 255.0;
}

//---------------------------------------------------------------------------
/*
 TriCairo_Color::TriCairo_Color(TColor t, double a)
 {
 R = (double)(t & 0x000000ff) / 255.0;
 G = (double)( (t & 0x0000ff00) >> 8) / 255.0;
 B = (double)( (t & 0x00ff0000) >> 16) / 255.0;
 A = a;
 }
 */

//---------------------------------------------------------------------------
void TriCairo_Color::Dispatch(double &r, double &g, double &b, double &a) {
  r = R;
  g = G;
  b = B;
  a = A;
}

//=============================================================================
//=============================================================================
//=============================================================================

TriCairo::TriCairo(unsigned int width, unsigned int height,
    TriCairo_CanvasType canvastype, String filename) :
    Width(width), Height(height) {
  Filename = filename;

  // Default constructor.
  CreateSurface(canvastype);
  cr = cairo_create(surface);

  // Clean up the surface.
  Clear(1.0, 1.0, 1.0);
  ColorD(0.0, 0.0, 0.0);

  //Bitmap = NULL;
}

//---------------------------------------------------------------------------
cairo_surface_t * TriCairo::CreateSurface(TriCairo_CanvasType canvastype) {
  CanvasType = canvastype;
  switch (CanvasType) {
    case ctBitmap:
      //surface = cairo_win32_surface_create_with_dib(CAIRO_FORMAT_ARGB32, Width,
      //    Height);
      break;
    case ctSurface:
      surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, Width, Height);
      break;
    case ctSVG:
      surface = cairo_svg_surface_create(Filename.c_str(), Width, Height);
      break;
    case ctPDF:
      surface = cairo_pdf_surface_create(Filename.c_str(), Width, Height);
      break;
    case ctPS:
      surface = cairo_ps_surface_create(Filename.c_str(), Width, Height);
      break;

  };
  return surface;
}

//---------------------------------------------------------------------------
TriCairo::~TriCairo(void) {
  /*
   if(Bitmap)
   {
   delete Bitmap;
   Bitmap = NULL;
   }
   */
  // Default destructor.
  if (cr) {
    cairo_destroy(cr);
    cr = NULL;
  }
  if (surface) {
    cairo_surface_destroy(surface);
    surface = NULL;
  }
}

//---------------------------------------------------------------------------
void TriCairo::Text(double x, double y, String Text,
    TriCairo_HorizontalAlignment ha, TriCairo_VerticalAlignment va) {
  double xo = 0.0;
  double yo = 0.0;

  cairo_text_extents_t extents;
  cairo_text_extents(cr, Text.c_str(), &extents);

  switch (ha) {
    case haLeft:
      xo = 0;
      break;
    case haCenter:
      xo = -extents.width / 2.0;
      break;
    case haRight:
      xo = -extents.width;
      break;
  }

  switch (va) {
    case vaBottom:
      yo = 0;
      break;
    case vaMiddle:
      yo = extents.height / 2.0;
      break;
    case vaTop:
      yo = extents.height;
      break;
  }

  cairo_move_to(cr, x + xo, y + yo);
  cairo_show_text(cr, Text.c_str());
}

//---------------------------------------------------------------------------
void TriCairo::Font(String Name, double Size, TriCairo_FontStyle Style) {
  switch (Style) {
    case Taquart::fsNormal:
      cairo_select_font_face(cr, Name.c_str(), CAIRO_FONT_SLANT_NORMAL,
          CAIRO_FONT_WEIGHT_NORMAL);
      break;
    case Taquart::fsItalic:
      cairo_select_font_face(cr, Name.c_str(), CAIRO_FONT_SLANT_ITALIC,
          CAIRO_FONT_WEIGHT_NORMAL);
      break;
    case Taquart::fsOblique:
      cairo_select_font_face(cr, Name.c_str(), CAIRO_FONT_SLANT_OBLIQUE,
          CAIRO_FONT_WEIGHT_NORMAL);
      break;
    case Taquart::fsNormalBold:
      cairo_select_font_face(cr, Name.c_str(), CAIRO_FONT_SLANT_NORMAL,
          CAIRO_FONT_WEIGHT_BOLD);
      break;
    case Taquart::fsItalicBold:
      cairo_select_font_face(cr, Name.c_str(), CAIRO_FONT_SLANT_ITALIC,
          CAIRO_FONT_WEIGHT_BOLD);
      break;
    case Taquart::fsObliqueBold:
      cairo_select_font_face(cr, Name.c_str(), CAIRO_FONT_SLANT_OBLIQUE,
          CAIRO_FONT_WEIGHT_BOLD);
      break;
  }

  cairo_set_font_size(cr, Size);
}

//---------------------------------------------------------------------------
void TriCairo::ColorD(double r, double g, double b, double a) {
  cairo_set_source_rgba(cr, r, g, b, a);
}

//---------------------------------------------------------------------------
void TriCairo::ColorB(unsigned int r, unsigned int g, unsigned int b,
    unsigned int a) {
  cairo_set_source_rgba(cr, double(r) / 255.0f, double(g) / 255.0f,
      double(b) / 255.0f, double(a) / 255.0f);
}

//---------------------------------------------------------------------------
void TriCairo::Color(TriCairo_Color Color) {
  double r = 0.0, g = 0.0, b = 0.0, a = 0.0;
  Color.Dispatch(r, g, b, a);
  cairo_set_source_rgba(cr, r, g, b, a);
}

//---------------------------------------------------------------------------
void TriCairo::Clear(double r, double g, double b, double a) {
  cairo_save(cr);
  if (a == 1.0f) {
    // Fully opaque.
    cairo_set_source_rgb(cr, r, g, b);
  }
  else {
    // Transparent.
    cairo_set_source_rgba(cr, r, g, b, a);
    cairo_set_operator(cr, CAIRO_OPERATOR_SOURCE);
  }
  cairo_paint(cr);
  cairo_restore(cr);
}

//---------------------------------------------------------------------------
void TriCairo::Clear(TriCairo_Color Color) {
  double r = 0.0, g = 0.0, b = 0.0, a = 0.0;
  Color.Dispatch(r, g, b, a);
  Clear(r, g, b, a);
}

//---------------------------------------------------------------------------
void TriCairo::MoveTo(double x, double y) {
  cairo_move_to(cr, x, y);
}

//---------------------------------------------------------------------------
void TriCairo::LineTo(double x, double y) {
  cairo_line_to(cr, x, y);
}

//---------------------------------------------------------------------------
void TriCairo::LineToRel(double x, double y) {
  cairo_rel_line_to(cr, x, y);
}

//---------------------------------------------------------------------------
void TriCairo::Arc(double x, double y, double r, double start, double end) {
  cairo_arc(cr, x, y, r, start, end);
}

void TriCairo::ClosePath(void) {
  cairo_close_path(cr);
}

//---------------------------------------------------------------------------
void TriCairo::Circle(double x, double y, double r) {
  Arc(x, y, r);
}

void TriCairo::Rectangle(double x, double y, double w, double h) {
  cairo_rectangle(cr, x, y, w, h);
}

//---------------------------------------------------------------------------
void TriCairo::LineWidth(double w) {
  cairo_set_line_width(cr, w);
}

//---------------------------------------------------------------------------
void TriCairo::Stroke(void) {
  cairo_stroke(cr);
}

//---------------------------------------------------------------------------
void TriCairo::StrokePreserve(void) {
  cairo_stroke_preserve(cr);
}

//---------------------------------------------------------------------------
void TriCairo::Fill(void) {
  cairo_fill(cr);
}

//---------------------------------------------------------------------------
void TriCairo::FillPreserve(void) {
  cairo_fill_preserve(cr);
}

//---------------------------------------------------------------------------
void TriCairo::LineCap(TriCairo_LineCap lc) {
  switch (lc) {
    case lcButt:
      cairo_set_line_cap(cr, CAIRO_LINE_CAP_BUTT); /* default */
      break;
    case lcRound:
      cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
      break;
    case lcSquare:
      cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
      break;
  }
}

//---------------------------------------------------------------------------
void TriCairo::LineJoin(TriCairo_LineJoin lj) {
  switch (lj) {
    case ljMiter:
      cairo_set_line_join(cr, CAIRO_LINE_JOIN_MITER); /* default */
      break;
    case ljBevel:
      cairo_set_line_join(cr, CAIRO_LINE_JOIN_BEVEL);
      break;
    case ljRound:
      cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);
      break;
  }
}

//---------------------------------------------------------------------------
void TriCairo::Save(String filename) {
  if (CanvasType == ctBitmap) {
    /*
     Graphics::TBitmap * Bitmap = new Graphics::TBitmap;
     Bitmap -> Width = Width;
     Bitmap -> Height = Height;
     BitBlt(Bitmap -> Canvas -> Handle, 0, 0, Width, Height,
     cairo_win32_surface_get_dc(surface), 0, 0, SRCCOPY);

     // Save result to output file.
     Bitmap -> SaveToFile(filename);
     delete Bitmap;
     */
  }
  else if (CanvasType == ctSurface) {
    cairo_surface_write_to_png(surface, filename.c_str());
  }
}

//=============================================================================
//=============================================================================
//=============================================================================

void Taquart::StrikeDipRake2MT(double S, double D, double R, double &M11,
    double &M22, double &M33, double &M12, double &M13, double &M23) {
  M11 = -1.0
      * (sin(D) * cos(R) * sin(2.0 * S)
          + sin(2.0 * D) * sin(R) * sin(S) * sin(S));
  M22 = (sin(D) * cos(R) * sin(2.0 * S)
      - sin(2.0 * D) * sin(R) * cos(S) * cos(S));
  M33 = -1.0 * (M11 + M22);
  M12 = (sin(D) * cos(R) * cos(2 * S)
      + 0.5 * sin(2.0 * D) * sin(R) * sin(2.0 * S)); /* Mxy */
  M13 = -1.0 * (cos(D) * cos(R) * cos(S) + cos(2 * D) * sin(R) * sin(S));
  M23 = -1.0 * (cos(D) * cos(R) * sin(S) - cos(2 * D) * sin(R) * cos(S));
}

//-----------------------------------------------------------------------------
double Taquart::zero_360(double str) {
  if (str >= 360.0)
    str -= 360.0;
  else if (str < 0.0)
    str += 360.0;
  return (str);
}

//-----------------------------------------------------------------------------
void Taquart::sincosd(double i, double *s, double *c) {
  *s = sin(i * DEG2RAD);
  *c = cos(i * DEG2RAD);
}

//-------------------------------------------------------------------------------------------------
// TRICAIRO
void Taquart::sincos(double i, double *s, double *c) {
  *s = sin(i);
  *c = cos(i);
}

//-------------------------------------------------------------------------------------------------
double Taquart::computed_dip1(nodal_plane NP1) {
  /* Genevieve Patau */
  /* Compute second nodal plane dip when are given strike, dip and rake for
   the first nodal plane with AKI & RICHARD's convention.
   Angles are in degrees. */

  double am = NP1.rake == 0.0 ? 1.0 : NP1.rake / fabs(NP1.rake);
  return acos(am * sind(NP1.rake) * sind(NP1.dip)) / DEG2RAD;

  // Last tested: $Id: utilmeca.c 12822 2014-01-31 23:39:56Z remko $ OK
}

/* GMT 5.1.2
 double computed_dip1 (struct nodal_plane NP1)
 {
 Compute second nodal plane dip when are given strike,
 dip and rake for the first nodal plane with AKI & RICHARD's
 convention.  Angles are in degrees.
 Genevieve Patau

 double am = (GMT_IS_ZERO (NP1.rake) ? 1.0 : NP1.rake / fabs (NP1.rake));
 double dip2;

 dip2 = acosd (am * sind (NP1.rake) * sind (NP1.dip));

 return (dip2);
 }
 */

//-------------------------------------------------------------------------------------------------
double Taquart::computed_rake1(struct nodal_plane NP1) {
  /* Compute rake in the second nodal plane when strike ,dip
   and rake are given for the first nodal plane with AKI &
   RICHARD's convention.
   Angles are in degrees. */
  /* Genevieve Patau */

  double sinrake2;
  double str2 = Taquart::computed_strike1(NP1);
  double dip2 = Taquart::computed_dip1(NP1);
  double am = NP1.rake == 0.0 ? 1.0 : NP1.rake / fabs(NP1.rake);
  double sd, cd, ss, cs;
  double sd2, cd2;
  Taquart::sincos(NP1.dip * DEG2RAD, &sd, &cd);
  Taquart::sincos(dip2 * DEG2RAD, &sd2, &cd2); /*  Onur 25Sep02,  */
  Taquart::sincos((NP1.str - str2) * DEG2RAD, &ss, &cs);

  if (fabs(dip2 - 90.0) < EPSIL)
    sinrake2 = am * cd;
  else
    sinrake2 = -am * sd * cs / cd2; /* ONUR: cd2 [cos(DIP2)] must be used not cd [cos(DIP1)] */
  // This part is different in GMT 5.1, as there cd is used instead of cd2

  return Taquart::datan2(sinrake2, -am * sd * ss);
}

/* GMT 5.1.2
 double computed_rake1 (struct nodal_plane NP1)
 {

 Compute rake in the second nodal plane when strike ,dip
 and rake are given for the first nodal plane with AKI &
 RICHARD's convention.  Angles are in degrees.

 Genevieve Patau


 double computed_strike1(), computed_dip1();
 double rake2, sinrake2;
 double str2 = computed_strike1(NP1);
 double dip2 = computed_dip1(NP1);
 double am = (GMT_IS_ZERO (NP1.rake) ? 1.0 : NP1.rake / fabs (NP1.rake));
 double sd, cd, ss, cs;
 sincosd (NP1.dip, &sd, &cd);
 sincosd (NP1.str - str2, &ss, &cs);

 if (fabs (dip2 - 90.0) < EPSIL)
 sinrake2 = am * cd;
 else
 sinrake2 = -am * sd * cs / cd;

 rake2 = d_atan2d (sinrake2, -am * sd * ss);

 return (rake2);
 }
 */

//-------------------------------------------------------------------------------------------------
double Taquart::computed_strike1(struct nodal_plane NP1) {
  /* Compute the strike of the decond nodal plane when are given strike, dip and rake for the first
   nodal plane with AKI & RICHARD's convention. Angles are in degrees. */
  /* Genevieve Patau */

  double str2;
  double cd1 = cosd(NP1.dip);
  double temp;
  double cp2, sp2;
  double am = NP1.rake == 0.0 ? 1.0 : NP1.rake / fabs(NP1.rake);
  double ss, cs, sr, cr;

  Taquart::sincos(NP1.rake * DEG2RAD, &sr, &cr);
  Taquart::sincos(NP1.str * DEG2RAD, &ss, &cs);
  if (cd1 < EPSIL && fabs(cr) < EPSIL) {
    str2 = NP1.str + 180.0;
  }
  else {
    temp = cr * cs;
    temp += sr * ss * cd1;
    sp2 = -am * temp;
    temp = ss * cr;
    temp -= sr * cs * cd1;
    cp2 = am * temp;
    str2 = datan2(sp2, cp2);
    str2 = Taquart::zero_360(str2);
  }
  return (str2);

}
/* GMT 5.1.2
 double computed_strike1 (struct nodal_plane NP1)
 {

 Compute the strike of the decond nodal plane when are given
 strike, dip and rake for the first nodal plane with AKI & RICHARD's
 convention.  Angles are in degrees.
 Genevieve Patau


 double str2, temp, cp2, sp2, ss, cs, sr, cr;
 double cd1 = cosd (NP1.dip);
 double am = (GMT_IS_ZERO (NP1.rake) ? 1. : NP1.rake /fabs (NP1.rake));

 sincosd (NP1.rake, &sr, &cr);
 sincosd (NP1.str, &ss, &cs);
 if (cd1 < EPSIL && fabs (cr) < EPSIL) {
 #if 0
 GMT_Report (API, GMT_MSG_DEBUG, "\nThe second plane is horizontal;");
 GMT_Report (API, GMT_MSG_DEBUG, "\nStrike is undetermined.");
 GMT_Report (API, GMT_MSG_DEBUG, "\nstr2 = NP1.str + 180. is taken to define");
 GMT_Report (API, GMT_MSG_DEBUG, "\nrake in the second plane.\n");
 #endif
 str2 = NP1.str + 180.0;
 }
 else {
 temp = cr * cs;
 temp += sr * ss * cd1;
 sp2 = -am * temp;
 temp = ss * cr;
 temp -= sr *  cs * cd1;
 cp2 = am * temp;
 str2 = d_atan2d(sp2, cp2);
 str2 = zero_360(str2);
 }
 return (str2);
 }
 */

//-------------------------------------------------------------------------------------------------
// ONUR
void Taquart::define_second_plane(struct nodal_plane NP1,
    struct nodal_plane *NP2) {
  /* Compute strike, dip, slip for the second nodal plane
   when are given strike, dip and rake for the first one. */
  /*  Genevieve Patau */

  NP2->str = Taquart::computed_strike1(NP1);
  NP2->dip = Taquart::computed_dip1(NP1);
  NP2->rake = Taquart::computed_rake1(NP1);
}

//-------------------------------------------------------------------------------------------------
// ONUR
double Taquart::null_axis_dip(double str1, double dip1, double str2,
    double dip2) {
  /* Compute null axis dip when strike and dip are given for each
   nodal plane. Angles are in degrees. */
  /* Genevieve Patau */

  double den = asin(sind(dip1) * sind(dip2) * sind(str1 - str2)) / DEG2RAD;
  if (den < 0.0)
    den = -den;
  return (den);
}

/* GMT 5.1.2
 double null_axis_dip (double str1, double dip1, double str2, double dip2)
 {
 compute null axis dip when strike and dip are given
 for each nodal plane.  Angles are in degrees.

 Genevieve Patau

 double den;

 den = asind (sind (dip1) * sind (dip2) * sind (str1 - str2));
 if (den < 0.) den = -den;
 return (den);
 }
 */

//-------------------------------------------------------------------------------------------------
double Taquart::null_axis_strike(double str1, double dip1, double str2,
    double dip2) {
  /* Compute null axis strike when strike and dip are given for each
   nodal plane. Angles are in degrees. */

  /* Genevieve Patau */
  double phn, cosphn, sinphn;
  double sd1, cd1, sd2, cd2, ss1, cs1, ss2, cs2;

  Taquart::sincos(dip1 * DEG2RAD, &sd1, &cd1);
  Taquart::sincos(dip2 * DEG2RAD, &sd2, &cd2);
  Taquart::sincos(str1 * DEG2RAD, &ss1, &cs1);
  Taquart::sincos(str2 * DEG2RAD, &ss2, &cs2);

  cosphn = sd1 * cs1 * cd2 - sd2 * cs2 * cd1;
  sinphn = sd1 * ss1 * cd2 - sd2 * ss2 * cd1;
  if (sind(str1 - str2) < 0.0) {
    cosphn = -cosphn;
    sinphn = -sinphn;
  }
  phn = datan2(sinphn, cosphn);
  if (phn < 0.0)
    phn += 360.0;
  return (phn);
}

/* GMT 5.1.2
 double null_axis_strike (double str1, double dip1, double str2, double dip2)
 {

 Compute null axis strike when strike and dip are given
 for each nodal plane.   Angles are in degrees.

 Genevieve Patau


 double phn, cosphn, sinphn, sd1, cd1, sd2, cd2, ss1, cs1, ss2, cs2;

 sincosd (dip1, &sd1, &cd1);
 sincosd (dip2, &sd2, &cd2);
 sincosd (str1, &ss1, &cs1);
 sincosd (str2, &ss2, &cs2);

 cosphn = sd1 * cs1 * cd2 - sd2 * cs2 * cd1;
 sinphn = sd1 * ss1 * cd2 - sd2 * ss2 * cd1;
 if (sind(str1 - str2) < 0.0) {
 cosphn = -cosphn;
 sinphn = -sinphn;
 }
 phn = d_atan2d(sinphn, cosphn);
 if (phn < 0.0) phn += 360.0;
 return (phn);
 }
 */

//-------------------------------------------------------------------------------------------------
double Taquart::computed_rake2(double str1, double dip1, double str2,
    double dip2, double fault) {
  /* Genevieve Patau */
  /* Compute rake in the second nodal plane when strikes and dips for the first
   and the second nodal plane are given with an additional variabled
   characterizing the fault: +1.0 for the inverse fault, -1.0 for the normal
   fault. Angles are in degrees. */

  double sinrake2;
  double sd, cd2, ss, cs;
  double cd;

  Taquart::sincos((str1 - str2) * DEG2RAD, &ss, &cs);
  sd = sind(dip1);
  cd = cosd(dip1); /* cos(dip1) not dip2   Onur Tan*/
  cd2 = cosd(dip2);

  if (fabs(dip2 - 90.0) < EPSIL)
    sinrake2 = fault * cd; // This part is different in GMT as cos(dip2) is used instead of cos(dip1)
  else
    sinrake2 = -fault * sd * cs / cd2; // This is the same as in GMT 5.1

  return Taquart::datan2(sinrake2, -fault * sd * ss);

}

/* GMT 5.1.2
 double computed_rake2 (double str1, double dip1, double str2, double dip2, double fault)
 {
 Compute rake in the second nodal plane when strike and dip
 for first and second nodal plane are given with a double
 characterizing the fault :
 +1. inverse fault
 -1. normal fault.
 Angles are in degrees.

 Genevieve Patau

 double rake2, sinrake2, sd, cd, ss, cs;

 sincosd (str1 - str2, &ss, &cs);

 sd = sind(dip1);        cd = cosd(dip2);
 if (fabs (dip2 - 90.0) < EPSIL)
 sinrake2 = fault * cd;
 else
 sinrake2 = -fault * sd * cs / cd;

 rake2 = d_atan2d (sinrake2, - fault * sd * ss);

 return (rake2);
 }
 */

//-------------------------------------------------------------------------------------------------
/* GMT 5.1.2
 void dc2axe (st_me meca, struct AXIS *T, struct AXIS *N, struct AXIS *P)
 {

 From FORTRAN routines of Anne Deschamps :
 compute azimuth and plungement of P-T axis
 from nodal plane strikes, dips and rakes.


 double cd1, sd1, cd2, sd2, cp1, sp1, cp2, sp2;
 double amz, amx, amy, dx, px, dy, py;

 cd1 = cosd (meca.NP1.dip) * M_SQRT2;
 sd1 = sind (meca.NP1.dip) * M_SQRT2;
 cd2 = cosd (meca.NP2.dip) * M_SQRT2;
 sd2 = sind (meca.NP2.dip) * M_SQRT2;
 cp1 = - cosd (meca.NP1.str) * sd1;
 sp1 = sind (meca.NP1.str) * sd1;
 cp2 = - cosd (meca.NP2.str) * sd2;
 sp2 = sind (meca.NP2.str) * sd2;

 amz = - (cd1 + cd2);
 amx = - (sp1 + sp2);
 amy = cp1 + cp2;
 dx = atan2d (hypot(amx, amy), amz) - 90.0;
 px = atan2d (amy, -amx);
 if (px < 0.0) px += 360.0;
 if (dx < EPSIL) {
 if (px > 90.0 && px < 180.0) px += 180.0;
 if (px >= 180.0 && px < 270.0) px -= 180.0;
 }

 amz = cd1 - cd2;
 amx = sp1 - sp2;
 amy = - cp1 + cp2;
 dy = atan2d (hypot(amx, amy), -fabs(amz)) - 90.0;
 py = atan2d (amy, -amx);
 if (amz > 0.0) py -= 180.0;
 if (py < 0.0) py += 360.0;
 if (dy < EPSIL) {
 if (py > 90.0 && py < 180.0) py += 180.0;
 if (py >= 180.0 && py < 270.0) py -= 180.0;
 }

 if (meca.NP1.rake > 0.0) {
 P->dip = dy; P->str = py;
 T->dip = dx; T->str = px;
 }
 else {
 P->dip = dx; P->str = px;
 T->dip = dy; T->str = py;
 }

 N->str = null_axis_strike (T->str, T->dip, P->str, P->dip);
 N->dip = null_axis_dip (T->str, T->dip, P->str, P->dip);
 }
 */

/*
 void Trinity::dc_to_axe(st_me meca,struct AXIS *T,struct AXIS *N,struct AXIS *P)
 {
 // From FORTRAN routines of Anne Deschamps : compute azimuth and
 //   plungement of P-T axis from nodal plane strikes, dips and rakes.

 int im;
 int pure_strike_slip = 0;
 double cd1, sd1, cd2, sd2;
 double cp1, sp1, cp2, sp2;
 double amz, amx, amy, dx, px, dy, py;

 if(fabs(sind(meca.NP1.rake)) > EPSIL) im = (int) (meca.NP1.rake / fabs(meca.NP1.rake));
 else if(fabs(sind(meca.NP2.rake)) > EPSIL) im = (int) (meca.NP2.rake / fabs(meca.NP2.rake));
 else pure_strike_slip = 1;

 if(pure_strike_slip)
 {
 if(cosd(meca.NP1.rake) < 0.)
 {
 P->str = zero_360(meca.NP1.str + 45.);
 T->str = zero_360(meca.NP1.str - 45.);
 }
 else
 {
 P->str = zero_360(meca.NP1.str - 45.);
 T ->str = zero_360(meca.NP1.str + 45.);
 }
 P->dip = 0.;
 T->dip = 0.;
 }
 else
 {
 cd1 = cosd(meca.NP1.dip) *  M_SQRT2;
 sd1 = sind(meca.NP1.dip) *  M_SQRT2;
 cd2 = cosd(meca.NP2.dip) *  M_SQRT2;
 sd2 = sind(meca.NP2.dip) *  M_SQRT2;
 cp1 = - cosd(meca.NP1.str) * sd1;
 sp1 = sind(meca.NP1.str) * sd1;
 cp2 = - cosd(meca.NP2.str) * sd2;
 sp2 = sind(meca.NP2.str) * sd2;

 amz = - (cd1 + cd2);
 amx = - (sp1 + sp2);
 amy = cp1 + cp2;
 dx = atan2(sqrt(amx * amx + amy * amy), amz) - M_PI_2;
 px = atan2(amy, - amx);
 if(px < 0.)
 px += TWO_PI;

 amz = cd1 - cd2;
 amx = sp1 - sp2;
 amy = - cp1 + cp2;
 dy = atan2(sqrt(amx * amx + amy * amy), - fabs(amz)) - M_PI_2;
 py = atan2(amy, - amx);
 if(amz > 0.)
 py -= M_PI;

 if(py < 0.)
 py += TWO_PI;

 if(im == 1)
 {
 P->dip = dy;
 P->str = py;
 T->dip = dx;
 T->str = px;
 }
 else
 {
 P->dip = dx;
 P->str = px;
 T->dip = dy;
 T->str = py;
 }
 }

 T->str /= D2R;
 T->dip /= D2R;
 P->str /= D2R;
 P->dip /= D2R;

 N->str =  null_axis_strike(T->str, T->dip, P->str, P->dip);
 N->dip =  null_axis_dip(T->str, T->dip, P->str, P->dip);
 }


 //-------------------------------------------------------------------------------------------------
 */
void Taquart::axe2dc(AXIS T, AXIS P, nodal_plane *NP1, nodal_plane *NP2) {
  double pp, dp, pt, dt;
  double p1, d1, p2, d2;
  double PII = M_PI * 2.;
  double cdp, sdp, cdt, sdt;
  double cpt, spt, cpp, spp;
  double amz, amy, amx;
  double im;

  pp = P.str * DEG2RAD;
  dp = P.dip * DEG2RAD;
  pt = T.str * DEG2RAD;
  dt = T.dip * DEG2RAD;

  sincos(dp, &sdp, &cdp);
  sincos(dt, &sdt, &cdt);
  sincos(pt, &spt, &cpt);
  sincos(pp, &spp, &cpp);

  cpt *= cdt;
  spt *= cdt;
  cpp *= cdp;
  spp *= cdp;

  amz = sdt + sdp;
  amx = spt + spp;
  amy = cpt + cpp;
  d1 = atan2(sqrt(amx * amx + amy * amy), amz);
  p1 = atan2(amy, -amx);
  if (d1 > M_PI_2) {
    d1 = M_PI - d1;
    p1 += M_PI;
    if (p1 > PII)
      p1 -= PII;
  }
  if (p1 < 0.)
    p1 += PII;

  amz = sdt - sdp;
  amx = spt - spp;
  amy = cpt - cpp;
  d2 = atan2(sqrt(amx * amx + amy * amy), amz);
  p2 = atan2(amy, -amx);
  if (d2 > M_PI_2) {
    d2 = M_PI - d2;
    p2 += M_PI;
    if (p2 > PII)
      p2 -= PII;
  }
  if (p2 < 0.)
    p2 += PII;

  NP1->dip = d1 / DEG2RAD;
  NP1->str = p1 / DEG2RAD;
  NP2->dip = d2 / DEG2RAD;
  NP2->str = p2 / DEG2RAD;

  im = 1;
  if (dp > dt)
    im = -1;
  NP1->rake = computed_rake2(NP2->str, NP2->dip, NP1->str, NP1->dip, im);
  NP2->rake = computed_rake2(NP1->str, NP1->dip, NP2->str, NP2->dip, im);
}

//=============================================================================
//=============================================================================
//=============================================================================

#define RAD2DEG (180.0/M_PI)
#define EPSIL 0.0001
#define NP (4) // ????
#define NR_END 1
#define FREE_ARG char*
#define NR_END 1
#define FREE_ARG char*

#undef DEBUG_MODE
#ifdef DEBUG_MODE
#include <fstream>
#endif
//---------------------------------------------------------------------------
Taquart::TriCairo_Meca::TriCairo_Meca(unsigned int width, unsigned int height,
    TriCairo_CanvasType canvastype, String filename) :
    TriCairo(width, height, canvastype, filename) {
  Margin = Width * 0.01;
  BWidth = Width - Margin * 2;
  BHeight = Height - Margin * 2;
  BRadius = BWidth > BHeight ? BHeight / 2 : BWidth / 2;
  BXo = Margin + BWidth / 2;
  BYo = Margin + BHeight / 2;
  BOutlineWidth = double(BRadius) / 100.0;
  BOutlineColor = TCColor(0.0, 0.0, 0.0);
  BPlusColor = TCColor(0.85, 0.85, 0.85);
  BMinusColor = TCColor(0.99, 0.99, 0.99);
  BTensorOutline = TCColor(0.5, 0.5, 0.5);
  BDCColor = TCColor(0.0, 0.0, 0.0);
  BDCWidth = double(BRadius) / 130.0;

  AxisFontFace = "Arial";
  AxisFontSize = double(BRadius) / 6.0;

  StationMarkerSize = double(BRadius) / 30.0;
  StationFontFace = "Arial";
  StationFontSize = double(BRadius) / 15.0;
  StationMarkerType = mtBWCircle;
  StationPlusColor = TCColor(1.0, 0.0, 0.0);
  StationMinusColor = TCColor(0.0, 0.0, 1.0);
  StationTextColor = TCColor(0.0, 0.0, 0.0);

  Projection = prSchmidt;
  Hemisphere = heLower;

  DrawCross = true;
  DrawAxis = true;
  DrawStations = true;
  DrawDC = true;
  DrawStationName = true;
}

//---------------------------------------------------------------------------
Taquart::TriCairo_Meca::~TriCairo_Meca(void) {
}

//---------------------------------------------------------------------------
void Taquart::TriCairo_Meca::ps_circle(double x, double y, double r,
    TCColor FillColor) {
  // Draw a circle and fill it with shade color.
  LineWidth(BOutlineWidth);
  Circle(x, y, r);
  Color(BOutlineColor);
  StrokePreserve();
  Color(FillColor);
  Fill();
}

//---------------------------------------------------------------------------
void Taquart::TriCairo_Meca::Polygon(double x[], double y[], int npoints,
    TCColor OutlineColor, bool fill, TCColor FillColor, double OutlineWidth) {
  //==== Project according to the coordination system.
  Project(x, y, npoints);
  //==== END: Project according to the coordination system.

  // Draw polygon and eventually fill it with brush color.
  //MoveTo(x[0], y[0]);
  MoveTo(x[0], Height - y[0]);
  for (int i = 1; i < npoints; i++) {
    //LineTo(x[i], y[i]);
    LineTo(x[i], Height - y[i]);
  }
  LineWidth(OutlineWidth);
  Color(OutlineColor);
  if (fill) {
    StrokePreserve();
    Color(FillColor);
    Fill();
  }
  else {
    Stroke();
  }
}

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
void Taquart::TriCairo_Meca::Axis(AXIS A, String AxisText) {
  double xp, yp;
  double radius;
  double spp, cpp;

  sincosd(A.str, &spp, &cpp);
  radius = sqrt(1.0 - sin(A.dip * DEG2RAD));
  if (radius >= 0.97)
    radius = 0.97;
  xp = radius * spp * BRadius + BXo;
  yp = Height - (radius * cpp * BRadius + BYo); // Flip upside down

  //==== Project according to the coordination system.
  Project(xp, yp);
  //==== END: Project according to the coordination system.

  //==== Draw P & T axis.
  ColorB(0, 0, 0);
  Font(AxisFontFace, AxisFontSize, fsNormal);
  Text(xp, yp, AxisText, haCenter, vaMiddle);
  Stroke();
}

//---------------------------------------------------------------------------
void Taquart::TriCairo_Meca::Station(double Azimuth, double Takeoff,
    double Disp, String Label, double &mx, double &my, double error) {
  double GA[4];
  if (Takeoff == 90.0f)
    Takeoff = 89.75f;
  GA[3] = cos(Takeoff * DEG2RAD);
  const double help = sqrt(1.0f - GA[3] * GA[3]);
  GA[1] = cos(Azimuth * DEG2RAD) * help;
  GA[2] = sin(Azimuth * DEG2RAD) * help;
  Station(GA, Disp, Label, mx, my);
}

//---------------------------------------------------------------------------
void Taquart::TriCairo_Meca::Station(double GA[], double Disp, String Label,
    double &mxx, double &myy, double error) {
  double alfa = 0.0;
  if (GA[3] < 0.0)
    alfa = (M_PI - acos(GA[3])) / 2.0;
  else if (GA[3] >= 0.0)
    alfa = acos(GA[3]) / 2.0;

  double pd = 0.0;
  switch (Projection) {
    case prWulff:
      pd = tan(alfa);
      break;
    case prSchmidt:
      pd = sin(alfa) * sqrt(2.0);
      break;
  }

  const double skal = sqrt(GA[1] * GA[1] + GA[2] * GA[2]);

  double x = 0.0, y = 0.0;
  if (GA[3] < 0.0) {
    x = GA[1] * pd / skal;
    y = GA[2] * pd / skal;
  }
  else if (GA[3] > 0.0) {
    x = -GA[1] * pd / skal;
    y = -GA[2] * pd / skal;
  }

  double my = BYo + BRadius * x;
  double mx = Height - (BXo + BRadius * y); // Flip upside-down

  //==== Project according to the coordination system.
  Project(mx, my);
  //==== END: Project according to the coordination system.

  mxx = mx;
  myy = my;

  //==== Draw error circle
  if (error > 0.0) {
    LineWidth(1.0);
    Circle(mx, my, StationMarkerSize + 1.0 * error * BRadius);
    Color(TCColor(1.0, 0.0, 1.0, 0.5));
    Stroke();
  }

  //==== Draw station and marker.
  DrawStationMarker(mx, my, Disp, Label);
}

//---------------------------------------------------------------------------
void Taquart::TriCairo_Meca::CenterCross(void) {
  LineWidth(1.0);
  Color(TCColor(1.0, 0.0, 0.0));
  MoveTo(BXo, BYo - BRadius / 20.0);
  LineTo(BXo, BYo + BRadius / 20.0);
  MoveTo(BXo - BRadius / 20.0, BYo);
  LineTo(BXo + BRadius / 20.0, BYo);
  Stroke();
}

//---------------------------------------------------------------------------
void Taquart::TriCairo_Meca::DrawStationMarker(double mx, double my,
    double Disp, String Label) {
  double w = 0.0;
  double s = StationMarkerSize;
  switch (StationMarkerType) {
    case mtPlusMinus:
      w = StationMarkerSize / 3.0;
      if (Disp > 0.0) {
        LineWidth(1.0);
        MoveTo(mx - w, my - s);
        LineTo(mx + w, my - s);
        LineTo(mx + w, my - w);
        LineTo(mx + s, my - w);
        LineTo(mx + s, my + w);
        LineTo(mx + w, my + w);
        LineTo(mx + w, my + s);
        LineTo(mx - w, my + s);
        LineTo(mx - w, my + w);
        LineTo(mx - s, my + w);
        LineTo(mx - s, my - w);
        LineTo(mx - w, my - w);
        ClosePath();
        Color(StationPlusColor);
        FillPreserve();
        Color(TCColor(0.0, 0.0, 0.0));
      }
      else {
        LineWidth(1.0);
        Rectangle(mx - StationMarkerSize, my - w, 2.0 * StationMarkerSize,
            2.0 * w);
        Color(StationMinusColor);
        FillPreserve();
        Color(TCColor(0.0, 0.0, 0.0));
      }
      break;

    case mtCircle:
      if (Disp > 0.0) {
        LineWidth(1.0);
        Circle(mx, my, StationMarkerSize);
        Color(StationPlusColor);
        FillPreserve();
        Color(TCColor(0.0, 0.0, 0.0));
      }
      else {
        LineWidth(1.0);
        Circle(mx, my, StationMarkerSize);
        Color(StationMinusColor);
        FillPreserve();
        Color(TCColor(0.0, 0.0, 0.0));
      }
      break;

      // Draw coloured squares.
    case mtSquare:
      if (Disp > 0.0) {
        LineWidth(1.0);
        Rectangle(mx - StationMarkerSize, my - StationMarkerSize,
            StationMarkerSize * 2.0, StationMarkerSize * 2.0);
        Color(StationPlusColor);
        FillPreserve();
        Color(TCColor(0.0, 0.0, 0.0));
      }
      else {
        LineWidth(1.0);
        Rectangle(mx - StationMarkerSize, my - StationMarkerSize,
            StationMarkerSize * 2.0, StationMarkerSize * 2.0);
        Color(StationMinusColor);
        FillPreserve();
        Color(TCColor(0.0, 0.0, 0.0));
      }
      break;

      // Plot Filled or empty circle (classical way of showing C and T parts)
    case mtBWCircle:
      if (Disp > 0.0) {
        LineWidth(1.0);
        Circle(mx, my, StationMarkerSize);
        Color(StationPlusColor);
        FillPreserve();
        Color(TCColor(0.0, 0.0, 0.0));
      }
      else {
        LineWidth(1.0);
        Circle(mx, my, StationMarkerSize);
        Color(StationMinusColor);
        FillPreserve();
        Color(TCColor(0.0, 0.0, 0.0));
      }
      break;
  }
  Stroke();

  if (DrawStationName) {
    Color(StationTextColor);
    Font(StationFontFace, StationFontSize, fsNormal);
    Text(mx, my + StationMarkerSize * 1.3, Label, haCenter, vaTop);
    Stroke();
  }
}

//---------------------------------------------------------------------------
void Taquart::TriCairo_Meca::DoubleCouple(double strike, double dip) {
  // Plot double-couple lines.
  /* Originally by Genevieve Patau */

  int i = -1;
  GMT_LONG npoints;

  double x[500], y[500];
  double str, radius;
  const double increment = 2.0;
  double si, co;

  str = strike;
  while (str <= strike + 180.0 + EPSIL) {
    i++;
    radius = proj_radius2(strike, dip, str) * BRadius;
    sincosd(str, &si, &co);
    x[i] = BXo + radius * si;
//    y[i] = Height - (BYo + radius * co);
    y[i] = (BYo + radius * co);
    str += increment;
  }
  npoints = i + 1;

  // Draw double couple
  Polygon(x, y, npoints, BDCColor, false, TCColor(0.0, 0.0, 0.0), BDCWidth);
}

//---------------------------------------------------------------------------
void Taquart::TriCairo_Meca::Tensor(AXIS T, AXIS N, AXIS P) {
  int plot_zerotrace = 0;

  double x0 = BXo;
  double y0 = BYo;
  double radius_size = BRadius;

  int d, b = 1, m;
  int i, ii, n = 0, j = 1, j2 = 0, j3 = 0;
  GMT_LONG npoints;
  TCColor rgb1, rgb2;
  int big_iso = 0;

  double a[3], p[3], v[3];
  double vi, iso, f;
  double fir, s2alphan, alphan;
  double cfi, sfi, can, san;
  double cpd, spd, cpb, spb, cpm, spm;
  double cad, sad, cab, sab, cam, sam;
  double xz, xn, xe;
  double az = 0., azp = 0., takeoff, r;
  double azi[3][2];
  double x[400], y[400], x2[400], y2[400], x3[400], y3[400];
  double xp1[800], yp1[800], xp2[800], yp2[800];
  double si, co;
  int jp_flag;
  int djp, mjp;

#ifdef DEBUG_MODE
  std::ofstream DebugFile( (Filename + Taquart::String(".log")).c_str());
  DebugFile << "File start" << std::endl;
#endif

  // Cleanup points.
  for (int i = 0; i < 400; i++) {
    x[i] = 0.0;
    y[i] = 0.0;
    x2[i] = 0.0;
    y2[i] = 0.0;
    x3[i] = 0.0;
    y3[i] = 0.0;
  }

  for (int i = 0; i < 800; i++) {
    xp1[i] = 0.0;
    yp1[i] = 0.0;
    xp2[i] = 0.0;
    yp2[i] = 0.0;
  }

  a[0] = T.str;
  a[1] = N.str;
  a[2] = P.str;
  p[0] = T.dip;
  p[1] = N.dip;
  p[2] = P.dip;
  v[0] = T.val;
  v[1] = N.val;
  v[2] = P.val;

  vi = (v[0] + v[1] + v[2]) / 3.0;
  for (i = 0; i <= 2; i++)
    v[i] = v[i] - vi;

  if (fabs(squared(v[0]) + squared(v[1]) + squared(v[2])) < EPSIL) {
#ifdef DEBUG_MODE
    DebugFile << "Pure explosion" << std::endl;
#endif
    /* pure implosion-explosion */
    if (vi > 0.0) {
      ps_circle(x0, y0, radius_size, BPlusColor);
    }
    if (vi < 0.0) {
      ps_circle(x0, y0, radius_size, BMinusColor);
    }
    return;
  }

  if (fabs(v[0]) >= fabs(v[2])) {
    d = 0;
    m = 2;
  }
  else {
    d = 2;
    m = 0;
  }

  if (plot_zerotrace)
    vi = 0.;

  f = -v[1] / v[d];
  iso = vi / v[d];
  jp_flag = 0; // Added 16.02.2015
  djp = -1; // Added 16.02.2015
  mjp = -1; // Added 16.02.2015

  /* Cliff Frohlich, Seismological Research letters,
   * Vol 7, Number 1, January-February, 1996
   * Unless the isotropic parameter lies in the range
   * between -1 and 1 - f there will be no nodes whatsoever */

  // No nodes!
  if (iso < -1.0) {
    ps_circle(x0, y0, radius_size, BMinusColor);
    return;
  }
  else if (iso > 1.0 - f) {
    ps_circle(x0, y0, radius_size, BPlusColor);
    return;
  }

  sincosd(p[d], &spd, &cpd);
  sincosd(p[b], &spb, &cpb);
  sincosd(p[m], &spm, &cpm);
  sincosd(a[d], &sad, &cad);
  sincosd(a[b], &sab, &cab);
  sincosd(a[m], &sam, &cam);

#ifdef DEBUG_MODE
  DebugFile << "Through angles" << std::endl;
#endif
  for (i = 0; i < 360; i++) {
    fir = (double) i * DEG2RAD;
    s2alphan = (2. + 2. * iso) / (3. + (1. - 2. * f) * cos(2. * fir));

    if (s2alphan > 1.) {
      big_iso++;

      // Modification according to J. Pesicek
      if (d == 0 && m == 2) {
        jp_flag = 1;
        djp = 2;
        mjp = 0;
      }

      if (d == 2 && m == 0) {
        jp_flag = 2;
        djp = 0;
        mjp = 2;
      }

      d = djp;
      m = mjp;

      f = -v[1] / v[d];
      iso = vi / v[d];
      s2alphan = (2. + 2. * iso) / (3. + (1. - 2. * f) * cos(2. * fir));
      sincosd(p[d], &spd, &cpd);
      sincosd(p[b], &spb, &cpb);
      sincosd(p[m], &spm, &cpm);
      sincosd(a[d], &sad, &cad);
      sincosd(a[b], &sab, &cab);
      sincosd(a[m], &sam, &cam);

      alphan = asin(sqrt(s2alphan));
      sincos(fir, &sfi, &cfi);
      sincos(alphan, &san, &can);

      xz = can * spd + san * sfi * spb + san * cfi * spm;
      xn = can * cpd * cad + san * sfi * cpb * cab + san * cfi * cpm * cam;
      xe = can * cpd * sad + san * sfi * cpb * sab + san * cfi * cpm * sam;

      if (fabs(xn) < EPSIL && fabs(xe) < EPSIL) {
        takeoff = 0.;
        az = 0.;
      }
      else {
        az = atan2(xe, xn);
        if (az < 0.)
          az += M_PI * 2.;
        takeoff = acos(xz / sqrt(xz * xz + xn * xn + xe * xe));
      }
      if (takeoff > M_PI_2) {
        takeoff = M_PI - takeoff;
        az += M_PI;
        if (az > M_PI * 2.)
          az -= M_PI * 2.;
      }

      //r = M_SQRT2 * sin(takeoff / 2.); /// PROJECTION????
      //---- Projection to Schmidt or Wullf net?
      switch (Projection) {
        case prWulff:
          r = tan(takeoff / 2.);
          break;
        case prSchmidt:
          r = M_SQRT2 * sin(takeoff / 2.); // Schmidt
          break;
        default:
          r = 0.0;
          break;
      }
      //---- END: Projection to Schmidt or Wullf net?

      sincos(az, &si, &co);
      if (i == 0) {
        azi[i][0] = az;
        x[i] = x0 + radius_size * r * si;
        y[i] = y0 + radius_size * r * co;
        azp = az;
      }
      else {
        if (fabs(fabs(az - azp) - M_PI) < DEG2RAD * 10.) {
          azi[n][1] = azp;
          azi[++n][0] = az;
        }
        if (fabs(fabs(az - azp) - M_PI * 2.) < DEG2RAD * 2.) {
          if (azp < az)
            azi[n][0] += M_PI * 2.;
          else
            azi[n][0] -= M_PI * 2.;
        }
        switch (n) {
          case 0:
            x[j] = x0 + radius_size * r * si;
            y[j] = y0 + radius_size * r * co;
            j++;
            break;
          case 1:
            x2[j2] = x0 + radius_size * r * si;
            y2[j2] = y0 + radius_size * r * co;
            j2++;
            break;
          case 2:
            x3[j3] = x0 + radius_size * r * si;
            y3[j3] = y0 + radius_size * r * co;
            j3++;
            break;
        }
        azp = az;
      }
    }

    // END: Modification according to J. Pesicek
    else {
      alphan = asin(sqrt(s2alphan));
      sincos(fir, &sfi, &cfi);
      sincos(alphan, &san, &can);

      xz = can * spd + san * sfi * spb + san * cfi * spm;
      xn = can * cpd * cad + san * sfi * cpb * cab + san * cfi * cpm * cam;
      xe = can * cpd * sad + san * sfi * cpb * sab + san * cfi * cpm * sam;

      if (fabs(xn) < EPSIL && fabs(xe) < EPSIL) {
        takeoff = 0.;
        az = 0.;
      }
      else {
        az = atan2(xe, xn);
        if (az < 0.)
          az += M_PI * 2.;
        takeoff = acos(xz / sqrt(xz * xz + xn * xn + xe * xe));
      }

      if (takeoff > M_PI_2) {
        takeoff = M_PI - takeoff;
        az += M_PI;
        if (az > M_PI * 2.)
          az -= M_PI * 2.;
      }

      //---- Projection to Schmidt or Wullf net?
      switch (Projection) {
        case prWulff:
          r = tan(takeoff / 2.);
          break;
        case prSchmidt:
          r = M_SQRT2 * sin(takeoff / 2.); // Schmidt
          break;
        default:
          r = 0.0;
          break;
      }
      //---- END: Projection to Schmidt or Wullf net?

      sincos(az, &si, &co);
      if (i == 0) {
        azi[i][0] = az;
        x[i] = x0 + radius_size * r * si;
        //y[i] = Height - (y0 + radius_size * r * co);
        y[i] = (y0 + radius_size * r * co);
        azp = az;
      }
      else {
        if (fabs(fabs(az - azp) - M_PI) < DEG2RAD * 10.) {
          azi[n][1] = azp;
          azi[++n][0] = az;
        }
        if (fabs(fabs(az - azp) - M_PI * 2.) < DEG2RAD * 2.) {
          if (azp < az)
            azi[n][0] += M_PI * 2.;
          else
            azi[n][0] -= M_PI * 2.;
        }
        switch (n) {
          case 0:
            x[j] = x0 + radius_size * r * si;
            //y[j] = Height - (y0 + radius_size * r * co);
            y[j] = (y0 + radius_size * r * co);
            j++;
            break;
          case 1:
            x2[j2] = x0 + radius_size * r * si;
            //y2[j2] = Height - (y0 + radius_size * r * co);
            y2[j2] = (y0 + radius_size * r * co);
            j2++;
            break;
          case 2:
            x3[j3] = x0 + radius_size * r * si;
            //y3[j3] = Height - (y0 + radius_size * r * co);
            y3[j3] = (y0 + radius_size * r * co);
            j3++;
            break;
        }
        azp = az;
      }
    }
  } // for loop

  azi[n][1] = az;

  if (v[1] < 0.) {
    rgb1 = BPlusColor;
    rgb2 = BMinusColor;
  }
  else {
    rgb1 = BMinusColor;
    rgb2 = BPlusColor;
  }

  if (!big_iso) {
    /* end patch to fix big_iso case plotting problems.  */
    ps_circle(x0, y0, radius_size, rgb2);
  }
  else if (jp_flag == 1) {
    ps_circle(x0, y0, radius_size, rgb1);
    //for (i=0;i<=2;i++) {rgb1[i] = e_rgb[i]; rgb2[i] = c_rgb[i];}
    rgb1 = BMinusColor;
    rgb2 = BPlusColor;
  }
  /* second case added. JP, DEC 2010 */
  else if (jp_flag == 2) {
    ps_circle(x0, y0, radius_size, rgb1);
    //for (i=0;i<=2;i++) {rgb2[i] = e_rgb[i]; rgb1[i] = c_rgb[i];}
    rgb2 = BMinusColor;
    rgb1 = BPlusColor;
  }
  //ps_circle(x0, y0, radius_size, rgb2);

  //return;
  //GMT_setfill (GMT, F1, false);
  switch (n) {
    case 0:
      for (i = 0; i < 360; i++) {
        xp1[i] = x[i];
        yp1[i] = y[i];
      }
      npoints = i;
      Polygon(xp1, yp1, npoints, BTensorOutline, true, rgb1);
#ifdef DEBUG_MODE
      DebugFile << "A Polygon(xp1, yp1, npoints, BTensorOutline, true, rgb1);" << std::endl;
      for(int kk=0; kk<npoints; kk++) {
        DebugFile << xp1[kk] << " " << yp1[kk] << std::endl;
      }
#endif

      break;

    case 1:
      for (i = 0; i < j; i++) {
        xp1[i] = x[i];
        yp1[i] = y[i];
      }

      if (azi[0][0] - azi[0][1] > M_PI)
        azi[0][0] -= M_PI * 2.;
      else if (azi[0][1] - azi[0][0] > M_PI)
        azi[0][0] += M_PI * 2.;

      if (azi[0][0] < azi[0][1]) {
        for (az = azi[0][1] - DEG2RAD; az > azi[0][0]; az -= DEG2RAD)
        {
          sincos(az, &si, &co);
          xp1[i] = x0 + radius_size * si;
          //yp1[i++] = Height - (y0 + radius_size * co);
          yp1[i++] = (y0 + radius_size * co);
        }
      }
      else {
        for (az = azi[0][1] + DEG2RAD; az < azi[0][0]; az += DEG2RAD)
        {
          sincos(az, &si, &co);
          xp1[i] = x0 + radius_size * si;
          //yp1[i++] = Height - (y0 + radius_size * co);
          yp1[i++] = (y0 + radius_size * co);
        }
      }
      npoints = i;
      Polygon(xp1, yp1, npoints, BTensorOutline, true, rgb1);
#ifdef DEBUG_MODE
      DebugFile << "B Polygon(xp1, yp1, npoints, BTensorOutline, true, rgb1);" << std::endl;
      for(int kk=0; kk<npoints; kk++) {
        DebugFile << xp1[kk] << " " << yp1[kk] << std::endl;
      }
#endif

      for (i = 0; i < j2; i++) {
        xp2[i] = x2[i];
        yp2[i] = y2[i];
      }

      if (azi[1][0] - azi[1][1] > M_PI)
        azi[1][0] -= M_PI * 2.;
      else if (azi[1][1] - azi[1][0] > M_PI)
        azi[1][0] += M_PI * 2.;

      if (azi[1][0] < azi[1][1]) {
        for (az = azi[1][1] - DEG2RAD; az > azi[1][0]; az -= DEG2RAD)
        {
          sincos(az, &si, &co);
          xp2[i] = x0 + radius_size * si;
          //yp2[i++] = Height - (y0 + radius_size * co);
          yp2[i++] = (y0 + radius_size * co);
        }
      }
      else {
        for (az = azi[1][1] + DEG2RAD; az < azi[1][0]; az += DEG2RAD)
        {
          sincos(az, &si, &co);
          xp2[i] = x0 + radius_size * si;
          //yp2[i++] = Height - (y0 + radius_size * co);
          yp2[i++] = (y0 + radius_size * co);
        }
      }
      npoints = i;
      Polygon(xp2, yp2, npoints, BTensorOutline, true, rgb1);
#ifdef DEBUG_MODE
      DebugFile << "C Polygon(xp2, yp2, npoints, BTensorOutline, true, rgb1);" << std::endl;
      for(int kk=0; kk<npoints; kk++) {
        DebugFile << xp2[kk] << " " << yp2[kk] << std::endl;
      }
#endif
      break;

    case 2:
      for (i = 0; i < j3; i++) {
        xp1[i] = x3[i];
        yp1[i] = y3[i];
      }
      for (ii = 0; ii < j; ii++) {
        xp1[i] = x[ii];
        yp1[i++] = y[ii];
      }

      /* Removed by Jeremy Pesicek
       if(big_iso)
       {
       for(ii=j2-1; ii>=0; ii--)
       {
       xp1[i] = x2[ii]; yp1[i++] = y2[ii];
       }
       npoints = i;
       Polygon(xp1, yp1, npoints, BTensorOutline, true, rgb1);
       break;
       }
       */

      if (azi[2][0] - azi[0][1] > M_PI)
        azi[2][0] -= M_PI * 2.;
      else if (azi[0][1] - azi[2][0] > M_PI)
        azi[2][0] += M_PI * 2.;

      if (azi[2][0] < azi[0][1])
        for (az = azi[0][1] - DEG2RAD; az > azi[2][0]; az -= DEG2RAD) {
          sincos(az, &si, &co);
          xp1[i] = x0 + radius_size * si;
          //yp1[i++] = Height - (y0 + radius_size * co);
          yp1[i++] = (y0 + radius_size * co);
        }
      else
        for (az = azi[0][1] + DEG2RAD; az < azi[2][0]; az += DEG2RAD) {
          sincos(az, &si, &co);
          xp1[i] = x0 + radius_size * si;
          //yp1[i++] = Height - (y0 + radius_size * co);
          yp1[i++] = (y0 + radius_size * co);
        }
      npoints = i;
      Polygon(xp1, yp1, npoints, BTensorOutline, true, rgb1);
#ifdef DEBUG_MODE
      DebugFile << "D Polygon(xp1, yp1, npoints, BTensorOutline, true, rgb1);" << std::endl;
      for(int kk=0; kk<npoints; kk++) {
        DebugFile << xp1[kk] << " " << yp1[kk] << std::endl;
      }
#endif
      for (i = 0; i < j2; i++) {
        xp2[i] = x2[i];
        yp2[i] = y2[i];
      }

      if (azi[1][0] - azi[1][1] > M_PI)
        azi[1][0] -= M_PI * 2.;
      else if (azi[1][1] - azi[1][0] > M_PI)
        azi[1][0] += M_PI * 2.;

      if (azi[1][0] < azi[1][1])
        for (az = azi[1][1] - DEG2RAD; az > azi[1][0]; az -= DEG2RAD)
        {
          sincos(az, &si, &co);
          xp2[i] = x0 + radius_size * si;
          //yp2[i++] = Height - (y0 + radius_size * co);
          yp2[i++] = (y0 + radius_size * co);
        }
      else
        for (az = azi[1][1] + DEG2RAD; az < azi[1][0]; az += DEG2RAD)
        {
          sincos(az, &si, &co);
          xp2[i] = x0 + radius_size * si;
          //yp2[i++] = Height - (y0 + radius_size * co);
          yp2[i++] = (y0 + radius_size * co);
        }
      npoints = i;
      Polygon(xp2, yp2, npoints, BTensorOutline, true, rgb1);
#ifdef DEBUG_MODE
      DebugFile << "E Polygon(xp2, yp2, npoints, BTensorOutline, true, rgb1);" << std::endl;
      for(int kk=0; kk<npoints; kk++) {
        DebugFile << xp2[kk] << " " << yp2[kk] << std::endl;
      }
#endif
      break;
  }
}

//---------------------------------------------------------------------------
void Taquart::TriCairo_Meca::Project(double &X, double &Y) {
  if (Hemisphere == heUpper) {
    X = 2.0 * BXo - X;
    Y = 2.0 * BYo - Y;
  }
}

//---------------------------------------------------------------------------
void Taquart::TriCairo_Meca::Project(double * X, double * Y,
    unsigned int npoints) {
  // Project according to coordinate system.
  if (Hemisphere == heUpper) {
    for (unsigned int i = 0; i < npoints; i++) {
      X[i] = 2.0 * BXo - X[i];
      Y[i] = 2.0 * BYo - Y[i];
    }
  }
}

//---------------------------------------------------------------------------
double Taquart::TriCairo_Meca::proj_radius2(double str1, double dip1,
    double str) {
  /*
   Compute the vector radius for a given strike,
   equal area projection, inferior sphere.
   Strike and dip of the plane are given.
   printf("\nVertical plane : strike is constant.");
   printf("\nFor ps_mechanism r == 1 for str = str1");
   printf("\n            else r == 0. is used.");
   */
  /* Genevieve Patau */

  double dip, r;

  if (fabs(dip1 - 90.) < EPSIL) {
    r = (fabs(str - str1) < EPSIL || fabs(str - str1 - 180.0) < EPSIL) ?
        1.0 : 0.0;
  }
  else {
    dip = atan(tan(dip1 * DEG2RAD) * sin((str - str1) * DEG2RAD));
    switch (Projection) {
      case prWulff:
        r = tan(M_PI_4 - dip / 2.0);
        break;
      case prSchmidt:
        r = sqrt(2.0) * sin(M_PI_4 - dip / 2.0);
        break;
      default:
        r = 0;
        break;
    }

  }
  return (r);
}

//---------------------------------------------------------------------------
void Taquart::TriCairo_Meca::axe2dc(AXIS T, AXIS P, nodal_plane *NP1,
    nodal_plane *NP2) {
  double pp, dp, pt, dt;
  double p1, d1, p2, d2;
  double PII = M_PI * 2.;
  double cdp, sdp, cdt, sdt;
  double cpt, spt, cpp, spp;
  double amz, amy, amx;
  double im;

  pp = P.str * DEG2RAD;
  dp = P.dip * DEG2RAD;
  pt = T.str * DEG2RAD;
  dt = T.dip * DEG2RAD;

  sincos(dp, &sdp, &cdp);
  sincos(dt, &sdt, &cdt);
  sincos(pt, &spt, &cpt);
  sincos(pp, &spp, &cpp);

  cpt *= cdt;
  spt *= cdt;
  cpp *= cdp;
  spp *= cdp;

  amz = sdt + sdp;
  amx = spt + spp;
  amy = cpt + cpp;
  d1 = atan2(sqrt(amx * amx + amy * amy), amz);
  p1 = atan2(amy, -amx);
  if (d1 > M_PI_2) {
    d1 = M_PI - d1;
    p1 += M_PI;
    if (p1 > PII)
      p1 -= PII;
  }
  if (p1 < 0.)
    p1 += PII;

  amz = sdt - sdp;
  amx = spt - spp;
  amy = cpt - cpp;
  d2 = atan2(sqrt(amx * amx + amy * amy), amz);
  p2 = atan2(amy, -amx);
  if (d2 > M_PI_2) {
    d2 = M_PI - d2;
    p2 += M_PI;
    if (p2 > PII)
      p2 -= PII;
  }
  if (p2 < 0.)
    p2 += PII;

  NP1->dip = d1 / DEG2RAD;
  NP1->str = p1 / DEG2RAD;
  NP2->dip = d2 / DEG2RAD;
  NP2->str = p2 / DEG2RAD;

  im = 1;
  if (dp > dt)
    im = -1;
  NP1->rake = computed_rake2(NP2->str, NP2->dip, NP1->str, NP1->dip, im);
  NP2->rake = computed_rake2(NP1->str, NP1->dip, NP2->str, NP2->dip, im);
}

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
double Taquart::TriCairo_Meca::squared(double v) {
  return pow(v, 2.0);
}

//---------------------------------------------------------------------------
#define MAX_SWEEPS 50

Taquart::GMT_LONG Taquart::TriCairo_Meca::GMT_jacobi(double *a, GMT_LONG *n,
    GMT_LONG *m, double *d, double *v, double *b, double *z, GMT_LONG *nrots) {
  /*
   *
   * Find eigenvalues & eigenvectors of a square symmetric matrix by Jacobi's
   * method.  Given A, find V and D such that A = V * D * V-transpose, with
   * V an orthogonal matrix and D a diagonal matrix.  The eigenvalues of A
   * are on diag(D), and the j-th column of V is the eigenvector corresponding
   * to the j-th diagonal element of D.  Returns 0 if OK, -1 if it fails to
   * converge in MAX_SWEEPS.
   *
   * a is sent as a square symmetric matrix, of size n, and row dimension m.
   * Only the diagonal and super-diagonal elements of a will be used, so the
   * sub-diagonal elements could be used to preserve a, or could have been
   * destroyed by an earlier attempt to form the Cholesky decomposition of a.
   * On return, the super-diagonal elements are destroyed.  The diagonal and
   * sub-diagonal elements are unchanged.
   * d is returned as an n-vector containing the eigenvalues of a, sorted
   * so that d[i] >= d[j] when i < j.  d = diag(D).
   * v is returned as an n by n matrix, V, with row dimension m, and the
   * columns of v are the eigenvectors corresponding to the values in d.
   * b is an n-vector of workspace, used to keep a copy of the diagonal
   * elements which is updated only after a full sweep.
   * z is an n-vector of workspace, used to accumulate the updates to
   * the diagonal values of a during each sweep.  This reduces round-
   * off problems.
   * nrots is the number of rotations performed.  Bounds on round-off error
   * can be estimated from this if desired.
   *
   * Numerical Details:
   * The basic algorithms is in many textbooks.  The idea is to make an
   * infinite series (which turns out to be at quadratically convergent)
   * of steps, in each of which A_new = P-transpose * A_old * P, where P is
   * a plane-rotation matrix in the p,q plane, through an angle chosen to
   * zero A_new(p,q) and A_new(q,p).  The sum of the diagonal elements
   * of A is unchanged by these operations, but the sum of squares of
   * diagonal elements of a is increased by 2 * |A_old(p,q)| at each step.
   * Although later steps make non-zero again the previously zeroed entries,
   * the sum of squares of diagonal elements increases with each rotation,
   * while the sum of squares of off-diagonals keeps decreasing, so that
   * eventually A_new is diagonal to machine precision.  This should
   * happen in a few (3 to 7) sweeps.
   *
   * If only the eigenvalues are wanted then there are faster methods, but
   * if all eigenvalues and eigenvectors are needed, then this method is
   * only somewhat slower than the fastest method (Householder tri-
   * diagonalization followed by symmetric QR iterations), and this method
   * is numerically extremely stable.
   *
   * C G J Jacobi ("Ueber ein leichtes Vefahren, die in der Theorie der
   * Saekularstoerungen vorkommenden Gelichungen numerisch aufzuloesen",
   * Crelle's Journal, v. 30, pp. 51--94, 1846) originally searched the
   * entire (half) matrix for the largest |A(p,q)| to select each step.
   * When the method was developed for machine computation (R T Gregory,
   * "Computing eigenvalues and eigenvectors of a symmetric matrix on
   * the ILLIAC", Math. Tab. and other Aids to Comp., v. 7, pp. 215--220,
   * 1953) it was done with a series of "sweeps" through the upper triangle,
   * visiting all p,q in turn.  Later, D A Pope and C Tompkins ("Maximizing
   * functions of rotations - experiments concerning speed of diagonalization
   * of symmetric matrices using Jacobi's method", J Assoc. Comput. Mach.
   * v. 4, pp. 459--466, 1957) introduced a variant that skips small
   * elements on the first few sweeps.  The algorithm here was given by
   * Heinz Rutishauser (1918--1970) and published in Numer. Math. v. 9,
   * pp 1--10, 1966, and in Linear Algebra (the Handbook for Automatic
   * Computation, v. II), by James Hardy Wilkinson and C. Reinsch (Springer-
   * Verlag, 1971).  It also appears in Numerical Recipes.
   *
   * This algorithm takes care to avoid round-off error in several ways.
   * First, although there are four values of theta in (-pi, pi] that
   * would zero A(p,q), there is only one with magnitude <= pi/4.
   * This one is used.  This is most stable, and also has the effect
   * that, if A_old(p,p) >= A_old(q,q) then A_new(p,p) > A_new(q,q).
   * Two copies of the diagonal elements are maintained in d[] and b[].
   * d[] is updated immediately in each rotation, and each new rotation
   * is computed based on d[], so that each rotation gets the benefit
   * of the previous ones.  However, z[] is also used to accumulate
   * the sum of all the changes in the diagonal elements during one sweep,
   * and z[] is used to update b[] after each sweep.  Then b is copied
   * to d.  In this way, at the end of each sweep, d is reset to avoid
   * accumulating round-off.
   *
   * This routine determines whether y is small compared to x by testing
   * if (fabs(y) + fabs(x) == fabs(x) ).  It is assumed that the
   * underflow which may occur here is nevertheless going to allow this
   * expression to be evaluated as TRUE or FALSE and execution to
   * continue.  If the run environment doesn't allow this, the routine
   * won't work properly.
   *
   * programmer:  W. H. F. Smith, 7 June, 1991.
   * Revised: PW: 12-MAR-1998 for GMT 3.1
   * Revision by WHF Smith, March 03, 2000, to speed up loop indexes.
   */
  GMT_LONG p, q, pp, pq, mp1, pm, qm, nsweeps, j, jm, i, k;
  double sum, threshold, g, h, t, theta, c, s, tau;

  /* Begin by initializing v, b, d, and z.  v = identity matrix,
   b = d = diag(a), and z = 0:  */

  for(int i=0; i<((*m) * (*n)); i++) {
    *(v+i) = 0.0;
  }

  for(int i=0; i<*n; i++) {
    *(z+i) = 0.0;
  }
  //memset((void *) v, 0, (size_t) ((*m) * (*n) * sizeof(double)));
  //memset((void *) z, 0, (size_t) ((*n) * sizeof(double)));

  mp1 = (*m) + 1;

  for (p = 0, pp = 0; p < (*n); p++, pp += mp1) {
    v[pp] = 1.0;
    b[p] = a[pp];
    d[p] = b[p];
  }

  /* End of initializations.  Set counters and begin:  */

  (*nrots) = 0;
  nsweeps = 0;

  while (nsweeps < MAX_SWEEPS) {

    /* Sum off-diagonal elements of upper triangle.  */
    sum = 0.0;
    for (q = 1, qm = (*m); q < (*n); q++, qm += (*m)) {
      for (p = 0, pq = qm; p < q; p++, pq++) {
        sum += fabs(a[pq]);
      }
    }

    /* Exit this loop (converged) when sum == 0.0  */
    if (sum == 0.0)
      break;

    /* If (nsweeps < 3) do only bigger elements;  else all  */
    threshold = (nsweeps < 3) ? 0.2 * sum / ((*n) * (*n)) : 0.0;

    /* Now sweep whole upper triangle doing Givens rotations:  */

    for (q = 1, qm = (*m); q < (*n); q++, qm += (*m)) {
      for (p = 0, pm = 0, pq = qm; p < q; p++, pm += (*m), pq++) {
        /* In 3/2000 I swapped order of these loops,
         to allow simple incrementing of pq  */

        if (a[pq] == 0.0)
          continue; /* New 3/2000  */

        g = 100.0 * fabs(a[pq]);

        /* After four sweeps, if g is small relative
         to a(p,p) and a(q,q), skip the
         rotation and set a(p,q) to zero.  */

        if ((nsweeps > 3) && ((fabs(d[p]) + g) == fabs(d[p]))
            && ((fabs(d[q]) + g) == fabs(d[q]))) {
          a[pq] = 0.0;
        }
        else if (fabs(a[pq]) > threshold) {

          h = d[q] - d[p];

          if (h == 0.0) {
            t = 1.0; /* This if block is new 3/2000  */
          }
          else if ((fabs(h) + g) == fabs(h)) {
            t = a[pq] / h;
          }
          else {
            theta = 0.5 * h / a[pq];
            t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
            if (theta < 0.0)
              t = -t;
          }

          c = 1.0 / sqrt(1.0 + t * t);
          s = t * c;
          tau = s / (1.0 + c);

          h = t * a[pq];
          z[p] -= h;
          z[q] += h;
          d[p] -= h;
          d[q] += h;
          a[pq] = 0.0;

          for (j = 0; j < p; j++) {
            g = a[j + pm];
            h = a[j + qm];
            a[j + pm] = g - s * (h + g * tau);
            a[j + qm] = h + s * (g - h * tau);
          }
          for (j = p + 1, jm = (*m) * (p + 1); j < q; j++, jm += (*m)) {
            g = a[p + jm];
            h = a[j + qm];
            a[p + jm] = g - s * (h + g * tau);
            a[j + qm] = h + s * (g - h * tau);
          }
          for (j = q + 1, jm = (*m) * (q + 1); j < (*n); j++, jm += (*m)) {
            g = a[p + jm];
            h = a[q + jm];
            a[p + jm] = g - s * (h + g * tau);
            a[q + jm] = h + s * (g - h * tau);
          }

          for (j = 0; j < (*n); j++) {
            g = v[j + pm];
            h = v[j + qm];
            v[j + pm] = g - s * (h + g * tau);
            v[j + qm] = h + s * (g - h * tau);
          }

          (*nrots)++;
        }
      }
    }

    /* End of one sweep of the upper triangle.  */

    nsweeps++;

    for (p = 0; p < (*n); p++) {
      b[p] += z[p]; /* Update the b copy of diagonal  */
      d[p] = b[p]; /* Replace d with b to reduce round-off error  */
      z[p] = 0.0; /* Clear z.  */
    }
  }

  /* Get here via break when converged, or when nsweeps == MAX_SWEEPS.
   Sort eigenvalues by insertion:  */

  for (i = 0; i < (*n) - 1; i++) {
    k = i;
    g = d[i];
    for (j = i + 1; j < (*n); j++) { /* Find max location  */
      if (d[j] >= g) {
        k = j;
        g = d[j];
      }
    }
    if (k != i) { /*  Need to swap value and vector  */
      d[k] = d[i];
      d[i] = g;
      p = i * (*m);
      q = k * (*m);
      for (j = 0; j < (*n); j++) {
        g = v[j + p];
        v[j + p] = v[j + q];
        v[j + q] = g;
      }
    }
  }

  /* Return 0 if converged; else print warning and return -1:  */

  return (0);
}

void Taquart::TriCairo_Meca::GMT_momten2axe(M_TENSOR mt, AXIS *T, AXIS *N,
    AXIS *P) {
  /* This version uses GMT_jacobi and does not suffer from the convert_matrix bug */
  GMT_LONG j, nrots;
  GMT_LONG np = 3;
  double *a, *d, *b, *z, *v;
  double az[3], pl[3];

  a = (double *) GMT_memory(NULL, (size_t) np * np, sizeof(double));
  d = (double *) GMT_memory(NULL, (size_t) np, sizeof(double));
  b = (double *) GMT_memory(NULL, (size_t) np, sizeof(double));
  z = (double *) GMT_memory(NULL, (size_t) np, sizeof(double));
  v = (double *) GMT_memory(NULL, (size_t) np * np, sizeof(double));

  a[0] = mt.f[0];
  a[1] = mt.f[3];
  a[2] = mt.f[4];
  a[3] = mt.f[3];
  a[4] = mt.f[1];
  a[5] = mt.f[5];
  a[6] = mt.f[4];
  a[7] = mt.f[5];
  a[8] = mt.f[2];

  if (GMT_jacobi(a, &np, &np, d, v, b, z, &nrots))
    throw;

  for (j = 0; j < np; j++) {
    pl[j] = asin(-v[j * np]);
    az[j] = atan2(v[j * np + 2], -v[j * np + 1]);
    if (pl[j] <= 0.) {
      pl[j] = -pl[j];
      az[j] += M_PI;
    }
    if (az[j] < 0)
      az[j] += 2.0 * M_PI;
    else if (az[j] > 2.0 * M_PI)
      az[j] -= 2.0 * M_PI;
    pl[j] *= (180.0 / M_PI);
    az[j] *= (180.0 / M_PI);
  }
  T->val = d[0];  //T->e = mt.expo;
  T->str = az[0];
  T->dip = pl[0];
  N->val = d[1];  //N->e = mt.expo;
  N->str = az[1];
  N->dip = pl[1];
  P->val = d[2];  //P->e = mt.expo;
  P->str = az[2];
  P->dip = pl[2];

  GMT_free((void *) a);
  GMT_free((void *) d);
  GMT_free((void *) b);
  GMT_free((void *) z);
  GMT_free((void *) v);
}

void * Taquart::TriCairo_Meca::GMT_memory(void *prev_addr, GMT_LONG nelem,
    size_t size) {
  return calloc((size_t) nelem, size);
}

void Taquart::TriCairo_Meca::GMT_free(void *addr) {
  free(addr);
}

//=============================================================================
//=============================================================================
//=============================================================================

//=============================================================================
//=============================================================================
//=============================================================================

//=============================================================================
//=============================================================================
//=============================================================================

