//-----------------------------------------------------------------------------
// Source: inputdata.h
// Module: focimt
// Class for storing input data for moment tensor inversion.
//
// Copyright (c) 2013-2015, Grzegorz Kwiatek.
//
// Permission is hereby granted, free of charge, to any person or organization
// obtaining a copy of the software and accompanying documentation covered by
// this license (the "Software") to use, reproduce, display, distribute,
// execute, and transmit the Software, and to prepare derivative works of the
// Software, and to permit third-parties to whom the Software is furnished to
// do so, all subject to the following:
//
// The copyright notices in the Software and this entire statement, including
// the above license grant, this restriction and the following disclaimer,
// must be included in all copies of the Software, in whole or in part, and
// all derivative works of the Software, unless such copies or derivative
// works are solely in the form of machine-executable object code generated by
// a source language processor.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
// SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
// FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//-----------------------------------------------------------------------------
#ifndef inputdataH
#define inputdataH
//---------------------------------------------------------------------------
#include "moment_tensor.h"

namespace Taquart {
  //! Input data vector for the seismic moment tensor calculation.
  /*! This structure wraps the input data vector for the seismic moment tensor
   *  calculation while using usmtcore module.
   *  \ingroup foci
   */
  class SMTInputLine {
    public:
      // Internal element not exported to XML file
      int Key; /*!< Station key-id (not exported!). This value
       is used to distinguish between different phases
       marked on the same channel. */

      // Elements exported to XLM file
      Taquart::String Name; /*!< Station name.*/
      unsigned int Id; /*!< Station id number.*/
      Taquart::String Component; /*!< Component.*/
      Taquart::String MarkerType; /*!< Type of the marker used.*/
      double Start; /*!< Start time [s].*/
      double End; /*!< End time [s].*/
      double Duration; /*!< Duration of signal [s].*/
      double Displacement; /*!< Displacement [m]. */
      double Incidence; /*!< Angle of incidence [deg]. */
      double Azimuth; /*!< Azimuth between station and source [deg]. */
      double TakeOff; /*!< Takeoff angle [deg]. */
      double Distance; /*!< Distance between station and source [m]. */
      double Density; /*!< Density [km/m**3]. */
      double Velocity; /*!< Average velocity [m/s]. */
      bool PickActive;
      bool ChannelActive;

      //! Export data in XML format.
      //XMLNode xmlExport(XMLExporter &Exporter);
    private:
    protected:
  };

  //! Input data wrapper for the seismic moment tensor solution calculation.
  /*! This class is used internally in %Foci to pass the input data for the
   *  seismic moment tensor inversion to the Foci::SMTThreadStruct class and
   *  usmt core module. 
   *  \ingroup foci
   */
  class SMTInputData {
    public:
      //! Default constructor.
      SMTInputData(void);

      //! Default destructor.
      ~SMTInputData(void);

      //! Recalculate statistics.
      bool Recalculate(void);

      //! Retrieve input data statistics.
      /*! \param OMeanDuration Reference to the mean rupture time.
       *  \param OStdDuration Reference to the standard deviation of the
       *  rupture time data.
       *  \param OMeanDisplacement Reference to the mean displacement.
       *  \param OStdDisplacement Reference to the standard deviation of the
       *  displacement data.
       */
      void GetStatistics(double &OMeanDuration, double &OStdDuration, double &OMeanDisplacement,
          double &OStdDisplacement);

      //! Find specified channel
      /*! \param ChannelName Channel name.
       *  \param OIndex Index of the input data specified by the ChannelName
       *  variable.
       *  \return \p true if channel has been found, \p false if has not been
       *  founded.
       */
      bool Find(Taquart::String &ChannelName, unsigned int &OIndex);
      bool Find(int Key, unsigned int &OIndex);

      //! Clear data wrapper.
      void Clear(void);

      //! Add input data.
      /*! \param InputLine Input data to add.
       *  \return Current number of items.
       */
      unsigned int Add(Taquart::SMTInputLine &InputLine);

      //! Set rupture time
      /*! \param ARuptureTime Rupture time [seconds]
       */
      void AddRuptureTime(double ARuptureTime);

      //! Get rupture time.
      /*! \return Rupture time [seconds]
       */
      double GetRuptureTime(void);

      //! Get input data structure.
      /*! \param Index Index of line to get data from.
       *  \param InputLine Reference to input data structure.
       */
      void Get(unsigned int Index, Taquart::SMTInputLine &InputLine) throw (Taquart::TriEOutOfRange);

      void Set(unsigned int Index, Taquart::SMTInputLine &InputLine) throw (Taquart::TriEOutOfRange);

      //! Get displacement (~seismic moment) value.
      /*! \param Index Index of line to get data from it.
       *  \return Displacement value.
       */
      double GetDisplacement(const unsigned int &Index) throw (Taquart::TriEOutOfRange);

      //! Remove input data line.
      /*! \param Index Index of line to remove.
       */
      void Remove(unsigned int Index) throw (Taquart::TriEOutOfRange);

      //! Return number of lines in input data table.
      /*! \return Number of lines of input tablel.
       */
      unsigned int Count(void);

      //! Export data in XML format.
      //XMLNode xmlExport(XMLExporter &Exporter);

      //! Copy constructor.
      /*! \param ASource Reference to the source class.
       */
      SMTInputData(const SMTInputData &ASource);

      //! Assignment operator.
      /*! \param ASource Reference to the source class.
       */
      SMTInputData& operator=(const SMTInputData &ASource);

      //! Assign values from another Foci::SMTInputData class.
      /*! \param ASource Reference to the source class.
       */
      void Assign(const SMTInputData &ASource);

      //! Calculate rupture time.
      /*! \param Result \a TRUE when rupture time was calculated properly.
       */
      double CountRuptureTime(bool &Result);
    private:
      int Key;
      double MeanDuration;
      double StdDuration;
      double MeanDisplacement;
      double StdDisplacement;
      double RuptureTime; /*!< Averaged rupture time. */
      std::vector<Taquart::SMTInputLine> InputData;
    protected:
  };
}

//---------------------------------------------------------------------------
#endif
