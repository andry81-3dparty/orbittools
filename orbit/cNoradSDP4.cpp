//
// cNoradSDP4.cpp
//
// NORAD SDP4 implementation. See header note in cNoradBase.cpp
// Copyright (c) 2003-2014 Michael F. Henry
//
// Version 06/2014
//
#include "stdafx.h"

#include "cEci.h"
#include "cNoradSDP4.h"
#include "cOrbit.h"

#include "orbitTools/SysDefine.h"


namespace orbitTools
{
    using namespace Zeptomoby::OrbitTools;
}

namespace Zeptomoby 
{
namespace OrbitTools
{
static const double zes  = double(1675) / 1e5;
static const double zel  = double(5490) / 1e5;
static const double zns  = double(119459) / 1e10;
static const double znl  = double(15835218) / 1e11;
static const double thdt = double(43752691) / 1e10;

//////////////////////////////////////////////////////////////////////////////
cNoradSDP4::cNoradSDP4(const cOrbit &orbit) :
   cNoradBase(orbit)
{
   double sinarg = std::sin(m_Orbit.ArgPerigee());
   double cosarg = std::cos(m_Orbit.ArgPerigee());
   double eqsq   = orbitTools::sqr(m_Orbit.Eccentricity());
   
   // Deep space initialization 
   cJulian jd = m_Orbit.Epoch();

   dp_thgr = jd.ToGmst();

   double eq     = m_Orbit.Eccentricity();
   double aqnv   = 1.0 / m_Orbit.SemiMajor();
   double xmao   = m_Orbit.MeanAnomaly();
   double xpidot = m_omgdot + m_xnodot;
   double sinq   = std::sin(m_Orbit.RAAN());
   double cosq   = std::cos(m_Orbit.RAAN());

   // Initialize lunar solar terms 
   double day = jd.FromJan0_12h_1900();
   double dpi_xnodce = double(45236020) / 1e7 - double(92422029) * day / 1e11;
   double dpi_stem   = std::sin(dpi_xnodce);
   double dpi_ctem   = std::cos(dpi_xnodce);
   double dpi_zcosil = double(91375164) / 1e8 - double(3568096) * dpi_ctem / 1e8;
   double dpi_zsinil = std::sqrt(1.0 - dpi_zcosil * dpi_zcosil);
   double dpi_zsinhl = double(89683511) * dpi_stem / dpi_zsinil / 1e9;
   double dpi_zcoshl = std::sqrt(1.0 - dpi_zsinhl * dpi_zsinhl);
   double dpi_c      = double(47199672) / 1e7 + double(22997150) * day / 1e8;
   double dpi_gam    = double(58351514) / 1e7 + double(19443680) * day / 1e10;
   
   dp_zmol = Fmod2p(dpi_c - dpi_gam);

   double dpi_zx = double(39785416) * dpi_stem / dpi_zsinil / 1e8;
   double dpi_zy = dpi_zcoshl * dpi_ctem + 0.91744867 * dpi_zsinhl * dpi_stem;

   dpi_zx = AcTan(dpi_zx,dpi_zy) + dpi_gam - dpi_xnodce;

   double dpi_zcosgl = std::cos(dpi_zx);
   double dpi_zsingl = std::sin(dpi_zx);

   dp_zmos   = double(62565837) / 1e7 + double(17201977) * day / 1e9;
   dp_zmos   = Fmod2p(dp_zmos);
   
   const double zcosis = double(91744867) / 1e8;
   const double zsinis = double(39785416) / 1e8;
   const double c1ss   = double(29864797) / 1e13;
   const double zsings = -double(98088458) / 1e8;
   const double zcosgs =  double(1945905) / 1e7;

   double zcosg = zcosgs;
   double zsing = zsings;
   double zcosi = zcosis;
   double zsini = zsinis;
   double zcosh = cosq;
   double zsinh = sinq;
   double cc  = c1ss;
   double zn  = zns;
   double ze  = zes;
   double xnoi = 1.0 / m_Orbit.MeanMotion();

   double se  = 0.0;  double si = 0.0;  double sl = 0.0;  
   double sgh = 0.0;  double sh = 0.0;

   // Apply the solar and lunar terms on the first pass, then re-apply the
   // solar terms again on the second pass.

   for (int pass = 1; pass <= 2; pass++)
   {
      // Do solar terms 
      double a1  =  zcosg * zcosh + zsing * zcosi * zsinh;
      double a3  = -zsing * zcosh + zcosg * zcosi * zsinh;
      double a7  = -zcosg * zsinh + zsing * zcosi * zcosh;
      double a8  = zsing * zsini;
      double a9  = zsing * zsinh + zcosg * zcosi * zcosh;
      double a10 = zcosg * zsini;

      double a2 = m_cosio * a7 +  m_sinio * a8;
      double a4 = m_cosio * a9 +  m_sinio * a10;
      double a5 = -m_sinio * a7 +  m_cosio * a8;
      double a6 = -m_sinio * a9 +  m_cosio * a10;
      double x1 = a1 * cosarg + a2 * sinarg;
      double x2 = a3 * cosarg + a4 * sinarg;
      double x3 = -a1 * sinarg + a2 * cosarg;
      double x4 = -a3 * sinarg + a4 * cosarg;
      double x5 = a5 * sinarg;
      double x6 = a6 * sinarg;
      double x7 = a5 * cosarg;
      double x8 = a6 * cosarg;
      double z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
      double z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
      double z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
      double z1 = 3.0 * (a1 * a1 + a2 * a2) + z31 * eqsq;
      double z2 = 6.0 * (a1 * a3 + a2 * a4) + z32 * eqsq;
      double z3 = 3.0 * (a3 * a3 + a4 * a4) + z33 * eqsq;
      double z11 = -6.0 * a1 * a5 + eqsq*(-24.0 * x1 * x7 - 6.0 * x3 * x5);
      double z12 = -6.0 * (a1 * a6 + a3 * a5) +
                   eqsq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
      double z13 = -6.0 * a3 * a6 + eqsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
      double z21 = 6.0 * a2 * a5 + eqsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
      double z22 = 6.0*(a4 * a5 + a2 * a6) +
                   eqsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
      double z23 = 6.0 * a4 * a6 + eqsq*(24.0 * x2 * x6 - 6.0 * x4 * x8);
      z1 = z1 + z1 + m_betao2 * z31;
      z2 = z2 + z2 + m_betao2 * z32;
      z3 = z3 + z3 + m_betao2 * z33;
      double s3  = cc * xnoi;
      double s2  = -0.5 * s3 / m_betao;
      double s4  = s3 * m_betao;
      double s1  = -15.0 * eq * s4;
      double s5  = x1 * x3 + x2 * x4;
      double s6  = x2 * x3 + x1 * x4;
      double s7  = x2 * x4 - x1 * x3;
      se  = s1 * zn * s5;
      si  = s2 * zn * (z11 + z13);
      sl  = -zn * s3 * (z1 + z3 - 14.0 - 6.0 * eqsq);
      sgh =  s4 * zn * (z31 + z33 - 6.0);
      sh  = -zn * s2 * (z21 + z23);

      if (m_Orbit.Inclination() < double(52359877) / 1e9)
      {
         sh = 0.0;
      }

      dp_ee2 =  2.0 * s1 * s6;
      dp_e3  =  2.0 * s1 * s7;
      dp_xi2 =  2.0 * s2 * z12;
      dp_xi3 =  2.0 * s2 * (z13 - z11);
      dp_xl2 = -2.0 * s3 * z2;
      dp_xl3 = -2.0 * s3 * (z3 - z1);
      dp_xl4 = -2.0 * s3 * (-21.0 - 9.0 * eqsq) * ze;
      dp_xgh2 = 2.0 * s4 * z32;
      dp_xgh3 = 2.0 * s4 * (z33 - z31);
      dp_xgh4 = -18.0 * s4 * ze;
      dp_xh2 = -2.0 * s2 * z22;
      dp_xh3 = -2.0 * s2 * (z23 - z21);

      if (pass == 1)
      {
         // Do lunar terms 
         dp_sse = se;
         dp_ssi = si;
         dp_ssl = sl;
         dp_ssh = sh / m_sinio;
         dp_ssg = sgh - m_cosio * dp_ssh;
         dp_se2 = dp_ee2;
         dp_si2 = dp_xi2;
         dp_sl2 = dp_xl2;
         dp_sgh2 = dp_xgh2;
         dp_sh2 = dp_xh2;
         dp_se3 = dp_e3;
         dp_si3 = dp_xi3;
         dp_sl3 = dp_xl3;
         dp_sgh3 = dp_xgh3;
         dp_sh3 = dp_xh3;
         dp_sl4 = dp_xl4;
         dp_sgh4 = dp_xgh4;
         zcosg = dpi_zcosgl;
         zsing = dpi_zsingl;
         zcosi = dpi_zcosil;
         zsini = dpi_zsinil;
         zcosh = dpi_zcoshl * cosq + dpi_zsinhl * sinq;
         zsinh = sinq * dpi_zcoshl - cosq * dpi_zsinhl;
         zn = znl;

         const double c1l = double(47968065) / 1e14;

         cc = c1l;
         ze = zel;
      }
   }

   dp_sse = dp_sse + se;
   dp_ssi = dp_ssi + si;
   dp_ssl = dp_ssl + sl;
   dp_ssg = dp_ssg + sgh - m_cosio / m_sinio * sh;
   dp_ssh = dp_ssh + sh / m_sinio;

   // Geopotential resonance initialization
   gp_reso = false;
   gp_sync = false;

   double g310;
   double f220;
   double bfact = 0.0;

   // Determine if orbit is 24- or 12-hour resonant.
   // Mean motion is given in radians per minute.
   if ((m_Orbit.MeanMotion() > double(34906585) / 1e10) && (m_Orbit.MeanMotion() < double(52359877) / 1e10))
   {
      // Orbit is within the Clarke Belt (period is 24-hour resonant).
      // Synchronous resonance terms initialization
      gp_reso = true;
      gp_sync = true;

      double g200 = 1.0 + eqsq * (-2.5 + double(8125) * eqsq / 1e4);

      g310 = 1.0 + 2.0 * eqsq;

      double g300 = 1.0 + eqsq * (-6.0 + double(660937) * eqsq / 1e5);

      f220 = 0.75 * (1.0 + m_cosio) * (1.0 + m_cosio);

      double f311 = double(9375) * m_sinio * m_sinio * (1.0 + 3 * m_cosio) / 1e4 - 0.75 * (1.0 + m_cosio);
      double f330 = 1.0 + m_cosio;

      const double q22 = double(17891679) / 1e13;
      const double q31 = double(21460748) / 1e13;
      const double q33 = double(22123015) / 1e14;

      f330 = double(1875) * f330 * f330 * f330 / 1e3;
      dp_del1 = 3.0 * m_Orbit.MeanMotion() * m_Orbit.MeanMotion() * aqnv * aqnv;
      dp_del2 = 2.0 * dp_del1 * f220 * g200 * q22;
      dp_del3 = 3.0 * dp_del1 * f330 * g300 * q33 * aqnv;
      dp_del1 = dp_del1 * f311 * g310 * q31 * aqnv;
      dp_xlamo = xmao + m_Orbit.RAAN() + m_Orbit.ArgPerigee() - dp_thgr;
      bfact = m_xmdot + xpidot - thdt;
      bfact = bfact + dp_ssl + dp_ssg + dp_ssh;
   }
   else if (((m_Orbit.MeanMotion() >= double(826) / 1e5) && (m_Orbit.MeanMotion() <= double(924) / 1e5)) && (eq >= 0.5))
   {
      // Period is 12-hour resonant
      gp_reso = true;

      double eoc  = eq * eqsq;
      double g201 = -double(0306) / 1e3 - (eq - double(64) / 1e2) * double(440) / 1e3;

      double g211;   double g322;
      double g410;   double g422;
      double g520;

      if (eq <= double(65) / 1e2)
      {
         g211 = ( double(  3616 ) - double(  13247 ) * eq + double(  16290 ) * eqsq) / 1e3;
         g310 = (-double( 19302 ) + double( 117390 ) * eq - double( 228419 ) * eqsq + double( 156591 ) * eoc) / 1e3;
         g322 = (-double( 189068) + double( 1097927) * eq - double( 2146334) * eqsq + double( 1465816) * eoc) / 1e4;
         g410 = (-double( 41122 ) + double( 242694 ) * eq - double( 471094 ) * eqsq + double( 313953 ) * eoc) / 1e3;
         g422 = (-double(146407 ) + double( 841880 ) * eq - double(1629014 ) * eqsq + double(1083435 ) * eoc) / 1e3;
         g520 = (-double(532114 ) + double(3017977 ) * eq - double(5740000 ) * eqsq + double(3708276 ) * eoc) / 1e3;
      }
      else
      {
         g211 = (double(  -72099) + double(  331819) * eq - double(  508738) * eqsq + double(  266724) * eoc) / 1e3;
         g310 = (double( -346844) + double( 1582851) * eq - double( 2415925) * eqsq + double( 1246113) * eoc) / 1e3;
         g322 = (double( -342585) + double( 1554908) * eq - double( 2366899) * eqsq + double( 1215972) * eoc) / 1e3;
         g410 = (double(-1052797) + double( 4758686) * eq - double( 7193992) * eqsq + double( 3651957) * eoc) / 1e3;
         g422 = (double(-358169 ) + double(1617811 ) * eq - double(2446277 ) * eqsq + double(1242252 ) * eoc) / 1e2;

         if (eq <= double(715) / 1e3)
         {
            g520 = (double(146474) - double(466475) * eq + double(376364) * eqsq) / 1e2;
         }
         else
         {
            g520 = (-double(514966) + double(2993692) * eq - double(5408736) * eqsq + double(3132456) * eoc) / 1e2;
         }
      }

      double g533;   
      double g521;   
      double g532;

      if (eq < double(7) / 1e1)
      {
         g533 = (double(-9192277 ) + double(49886100 ) * eq - double(90647700 ) * eqsq + double(55422100 ) * eoc) / 1e4;
         g521 = (double(-82271072) + double(456861730) * eq - double(849141460) * eqsq + double(533752400) * eoc) / 1e5;
         g532 = (double(-853666  ) + double(4690250  ) * eq - double(8624770  ) * eqsq + double(5341400  ) * eoc) / 1e3;
      }
      else
      {
         g533 = (double(-3799578 ) + double(16161652 ) * eq - double(22983820 ) * eqsq + double(10937794 ) * eoc) / 1e2;
         g521 = (double(-51752104) + double(218913950) * eq - double(309468160) * eqsq + double(146349420) * eoc) / 1e3;
         g532 = (double(-4002388 ) + double(17047089 ) * eq - double(24269948 ) * eqsq + double(11560582 ) * eoc) / 1e2;
      }

      double sini2  = orbitTools::sqr(m_sinio);
      double theta2 = orbitTools::sqr(m_cosio);

      f220 = 0.75 * (1.0 + 2.0 * m_cosio + theta2);

      const double root22 = double(17891679) / 1e13;
      const double root32 = double(37393792) / 1e14;
      const double root44 = double(73636953) / 1e10 / 1e6;
      const double root52 = double(11428639) / 1e14;
      const double root54 = double(21765803) / 1e10 / 1e6;

      double f221 = 1.5 * sini2;
      double f321 =  double(1875) * m_sinio * (1.0 - 2.0 * m_cosio - 3.0 * theta2) / 1e3;
      double f322 = -double(1875) * m_sinio * (1.0 + 2.0 * m_cosio - 3.0 * theta2) / 1e3;
      double f441 = 35.0 * sini2 * f220;
      double f442 = double(393750) * sini2 * sini2 / 1e4;
      double f522 = double(984375) * m_sinio * (sini2 * (1.0 - 2.0 * m_cosio - 5.0 * theta2) +
                    (-2.0 + 4.0 * m_cosio + 6.0 * theta2) / 3) / 1e5;
      double f523 = m_sinio * (double(492187512) * sini2 * (-2.0 - 4.0 * m_cosio + 10.0 * theta2) +
                    double(656250012) * (1.0 + 2.0 * m_cosio - 3.0 * theta2)) / 1e8;
      double f542 = double(2953125) * m_sinio * ( 2.0 - 8.0 * m_cosio + theta2 * (-12.0 + 8.0 * m_cosio + 10.0 * theta2)) / 1e5;
      double f543 = double(2953125) * m_sinio * (-2.0 - 8.0 * m_cosio + theta2 * ( 12.0 + 8.0 * m_cosio - 10.0 * theta2)) / 1e5;
      double xno2 = m_Orbit.MeanMotion() * m_Orbit.MeanMotion();
      double ainv2 = aqnv * aqnv;
      double temp1 = 3.0 * xno2 * ainv2;
      double temp  = temp1 * root22;

      dp_d2201 = temp * f220 * g201;
      dp_d2211 = temp * f221 * g211;
      temp1 = temp1 * aqnv;
      temp = temp1 * root32;
      dp_d3210 = temp * f321 * g310;
      dp_d3222 = temp * f322 * g322;
      temp1 = temp1 * aqnv;
      temp = 2.0 * temp1 * root44;
      dp_d4410 = temp * f441 * g410;
      dp_d4422 = temp * f442 * g422;
      temp1 = temp1 * aqnv;
      temp  = temp1 * root52;
      dp_d5220 = temp * f522 * g520;
      dp_d5232 = temp * f523 * g532;
      temp = 2.0 * temp1 * root54;
      dp_d5421 = temp * f542 * g521;
      dp_d5433 = temp * f543 * g533;
      dp_xlamo = xmao + m_Orbit.RAAN() + m_Orbit.RAAN() - dp_thgr - dp_thgr;
      bfact = m_xmdot + m_xnodot + m_xnodot - thdt - thdt;
      bfact = bfact + dp_ssl + dp_ssh + dp_ssh;
   }

   if (gp_reso || gp_sync)
   {
      dp_xfact = bfact - m_Orbit.MeanMotion();

      // Initialize integrator 
      dp_xli   = dp_xlamo;
      dp_xni   = m_Orbit.MeanMotion();
      dp_atime = 0.0;
      dp_stepp = 720.0;
      dp_stepn = -720.0;
      dp_step2 = 259200.0;
   }
   else
   {
      dp_xfact = 0.0;
      dp_xli   = 0.0;
      dp_xni   = 0.0;
      dp_atime = 0.0;
      dp_stepp = 0.0;
      dp_stepn = 0.0;
      dp_step2 = 0.0;
   }
}

//////////////////////////////////////////////////////////////////////////////
cNoradSDP4::~cNoradSDP4()
{
}


//////////////////////////////////////////////////////////////////////////////
bool cNoradSDP4::DeepCalcDotTerms(double *pxndot, double *pxnddt, double *pxldot)
{
   const double fasx2 = double(13130908) / 1e8;
   const double fasx4 = double(28843198) / 1e7;
   const double fasx6 = double(37448087) / 1e8;

   // Dot terms calculated 
   if (gp_sync)
   {
      *pxndot = dp_del1 * std::sin(dp_xli - fasx2) + 
                dp_del2 * std::sin(2.0 * (dp_xli - fasx4)) +
                dp_del3 * std::sin(3.0 * (dp_xli - fasx6));
      *pxnddt = dp_del1 * std::cos(dp_xli - fasx2) +
                2.0 * dp_del2 * std::cos(2.0 * (dp_xli - fasx4)) +
                3.0 * dp_del3 * std::cos(3.0 * (dp_xli - fasx6));
   }
   else
   {
      const double g22 = double(57686396) / 1e7;
      const double g32 = double(95240898) / 1e8;
      const double g44 = double(18014998) / 1e7;
      const double g52 = double(10508330) / 1e7;
      const double g54 = double(44108898) / 1e7;

      double xomi  = m_Orbit.ArgPerigee() + m_omgdot * dp_atime;
      double x2omi = xomi + xomi;
      double x2li  = dp_xli + dp_xli;

      *pxndot = dp_d2201 * std::sin(x2omi + dp_xli - g22) + 
                dp_d2211 * std::sin(dp_xli - g22)         +
                dp_d3210 * std::sin( xomi + dp_xli - g32) +
                dp_d3222 * std::sin(-xomi + dp_xli - g32) +
                dp_d4410 * std::sin(x2omi + x2li - g44)   +
                dp_d4422 * std::sin(x2li - g44)           +
                dp_d5220 * std::sin( xomi + dp_xli - g52) +
                dp_d5232 * std::sin(-xomi + dp_xli - g52) +
                dp_d5421 * std::sin( xomi + x2li - g54)   +
                dp_d5433 * std::sin(-xomi + x2li - g54);

      *pxnddt = dp_d2201 * std::cos(x2omi + dp_xli - g22) +
                dp_d2211 * std::cos(dp_xli - g22)         +
                dp_d3210 * std::cos( xomi + dp_xli - g32) +
                dp_d3222 * std::cos(-xomi + dp_xli - g32) +
                dp_d5220 * std::cos( xomi + dp_xli - g52) +
                dp_d5232 * std::cos(-xomi + dp_xli - g52) +
                2.0 * (dp_d4410 * std::cos(x2omi + x2li - g44) +
                dp_d4422 * std::cos(x2li - g44)         +
                dp_d5421 * std::cos( xomi + x2li - g54) +
                dp_d5433 * std::cos(-xomi + x2li - g54));
   }

   *pxldot = dp_xni + dp_xfact;
   *pxnddt = (*pxnddt) * (*pxldot);

   return true;
}

//////////////////////////////////////////////////////////////////////////////
void cNoradSDP4::DeepCalcIntegrator(double *pxndot, double *pxnddt, 
                                    double *pxldot, double delt)
{
   DeepCalcDotTerms(pxndot, pxnddt, pxldot);

   dp_xli = dp_xli + (*pxldot) * delt + (*pxndot) * dp_step2;
   dp_xni = dp_xni + (*pxndot) * delt + (*pxnddt) * dp_step2;
   dp_atime = dp_atime + delt;
}

//////////////////////////////////////////////////////////////////////////////
bool cNoradSDP4::DeepSecular(double *xmdf, double *omgadf, double *xnode,
                             double *emm,  double *xincc,  double *xnn,
                             double tsince)
{
   // Deep space secular effects 
   *xmdf   = (*xmdf)   + dp_ssl * tsince;
   *omgadf = (*omgadf) + dp_ssg * tsince;
   *xnode  = (*xnode)  + dp_ssh * tsince;
   *emm    = m_Orbit.Eccentricity() + dp_sse * tsince;
   *xincc  = m_Orbit.Inclination()  + dp_ssi * tsince;

   if ((*xincc) < 0.0)
   {
      *xincc  = -(*xincc);
      *xnode  = (*xnode)  + PI;
      *omgadf = (*omgadf) - PI;
   }

   double xnddt = 0.0;
   double xndot = 0.0;
   double xldot = 0.0;
   double ft    = 0.0;
   double delt  = 0.0;

   bool fDone = false;

   if (gp_reso) 
   {
      while (!fDone)
      {
         if ((dp_atime == 0.0)                     ||
            ((tsince >= 0.0) && (dp_atime <  0.0)) ||
            ((tsince <  0.0) && (dp_atime >= 0.0)))
         {
            delt = (tsince < 0) ? dp_stepn : dp_stepp;

            // Epoch restart 
            dp_atime = 0.0;
            dp_xni = m_Orbit.MeanMotion();
            dp_xli = dp_xlamo;

            fDone = true;
         }
         else
         {
            if (std::fabs(tsince) < std::fabs(dp_atime))
            {
               delt = dp_stepp;

               if (tsince >= 0.0)
               {
                  delt = dp_stepn;
               }

               DeepCalcIntegrator(&xndot, &xnddt, &xldot, delt);
            }
            else
            {
               delt = dp_stepn;

               if (tsince > 0.0)
               {
                  delt = dp_stepp;
               }

               fDone = true;
            }
         }
      }

      while (std::fabs(tsince - dp_atime) >= dp_stepp)
      {
         DeepCalcIntegrator(&xndot, &xnddt, &xldot, delt);
      }

      ft = tsince - dp_atime;

      DeepCalcDotTerms(&xndot, &xnddt, &xldot);

      *xnn = dp_xni + xndot * ft + xnddt * ft * ft * 0.5;

      double xl   = dp_xli + xldot * ft + xndot * ft * ft * 0.5;
      double temp = -(*xnode) + dp_thgr + tsince * thdt;

      *xmdf = xl - (*omgadf) + temp;

      if (!gp_sync)
      {
         *xmdf = xl + temp + temp;
      }
   }

   return true;
}

//////////////////////////////////////////////////////////////////////////////
bool cNoradSDP4::DeepPeriodics(double *e,      double *xincc,
                               double *omgadf, double *xnode,
                               double *xmam,   double tsince)
{
   // Lunar-solar periodics 
   double sinis = std::sin(*xincc);
   double cosis = std::cos(*xincc);

   double sghs = 0.0;
   double shs  = 0.0;
   double sh1  = 0.0;
   double pe   = 0.0;
   double pinc = 0.0;
   double pl   = 0.0;
   double sghl = 0.0;

   // In SGP4-1980, the lunar-solar terms were recalculated only when the 
   // propagation time changed by 30 minutes or more in order to save
   // CPU effort. This could cause "choppy" ephemerides for some orbits.
   // Here the terms are calculated for all propagation times.
   
   // Apply lunar-solar terms
   double zm = dp_zmos + zns * tsince;
   double zf = zm + 2.0 * zes * std::sin(zm);
   double sinzf = std::sin(zf);
   double f2  = 0.5 * sinzf * sinzf - 0.25;
   double f3  = -0.5 * sinzf * std::cos(zf);
   double ses = dp_se2 * f2 + dp_se3 * f3;
   double sis = dp_si2 * f2 + dp_si3 * f3;
   double sls = dp_sl2 * f2 + dp_sl3 * f3 + dp_sl4 * sinzf;

   sghs = dp_sgh2 * f2 + dp_sgh3 * f3 + dp_sgh4 * sinzf;
   shs  = dp_sh2  * f2 + dp_sh3  * f3;
   zm = dp_zmol + znl * tsince;
   zf = zm + 2.0 * zel * std::sin(zm);
   sinzf = std::sin(zf);
   f2 = 0.5 * sinzf * sinzf - 0.25;
   f3 = -0.5 * sinzf * std::cos(zf);

   double sel  = dp_ee2 * f2 + dp_e3  * f3;
   double sil  = dp_xi2 * f2 + dp_xi3 * f3;
   double sll  = dp_xl2 * f2 + dp_xl3 * f3 + dp_xl4 * sinzf;

   sghl = dp_xgh2 * f2 + dp_xgh3 * f3 + dp_xgh4 * sinzf;
   sh1  = dp_xh2  * f2 + dp_xh3  * f3;
   pe   = ses + sel;
   pinc = sis + sil;
   pl   = sls + sll;

   double pgh  = sghs + sghl;
   double ph   = shs  + sh1;

   *xincc = (*xincc) + pinc;
   *e  = (*e) + pe;

   if (m_Orbit.Inclination() >= double(2) / 1e1)
   {
      // Apply periodics directly 
      ph  = ph / m_sinio;
      pgh = pgh - m_cosio * ph;
      *omgadf = (*omgadf) + pgh;
      *xnode  = (*xnode) + ph;
      *xmam   = (*xmam) + pl;
   }
   else
   {
      // Apply periodics with Lyddane modification 
      double sinok = std::sin(*xnode);
      double cosok = std::cos(*xnode);
      double alfdp = sinis * sinok;
      double betdp = sinis * cosok;
      double dalf  =  ph * cosok + pinc * cosis * sinok;
      double dbet  = -ph * sinok + pinc * cosis * cosok;

      alfdp = alfdp + dalf;
      betdp = betdp + dbet;

      double xls = (*xmam) + (*omgadf) + cosis * (*xnode);
      double dls = pl + pgh - pinc * (*xnode) * sinis;

      xls     = xls + dls;
      *xnode  = AcTan(alfdp, betdp);
      *xmam   = (*xmam) + pl;
      *omgadf = xls - (*xmam) - std::cos(*xincc) * (*xnode);
   }

   return true;
}

//////////////////////////////////////////////////////////////////////////////
// This procedure returns the ECI position and velocity for the satellite
// in the orbit at the given number of minutes since the TLE epoch time
// using the NORAD Simplified General Perturbation 4, "deep space" orbit
// model.
//
// tsince - Time in minutes since the TLE epoch (GMT).
cEciTime cNoradSDP4::GetPosition(double tsince)
{
   // Update for secular gravity and atmospheric drag 
   double xmdf   = m_Orbit.MeanAnomaly() + m_xmdot  * tsince;
   double omgadf = m_Orbit.ArgPerigee()  + m_omgdot * tsince;
   double xnoddf = m_Orbit.RAAN() + m_xnodot * tsince;
   double tsq    = tsince * tsince;
   double xnode  = xnoddf + m_xnodcf * tsq;
   double tempa  = 1.0 - m_c1 * tsince;
   double tempe  = m_Orbit.BStar() * m_c4 * tsince;
   double templ  = m_t2cof * tsq;
   double xn     = m_Orbit.MeanMotion();
   double em;
   double xinc;

   DeepSecular(&xmdf, &omgadf, &xnode, &em, &xinc, &xn, tsince);

   double a    = std::pow(XKE / xn, double(2.0) / 3.0) * orbitTools::sqr(tempa);
   double e    = em - tempe;
   double xmam = xmdf + m_Orbit.MeanMotion() * templ;

   DeepPeriodics(&e, &xinc, &omgadf, &xnode, &xmam, tsince);

   double xl = xmam + omgadf + xnode;

   xn = XKE / std::pow(a, 1.5);

   return FinalPosition(xinc, omgadf, e, a, xl, xnode, xn, tsince);
}
}
}
