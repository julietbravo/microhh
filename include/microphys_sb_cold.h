/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "constants.h"
#include "fast_math.h"
#include "thermo_moist_functions.h"

namespace Sb_cold
{
    /*
       Kernels ported from ICON.
       All kernels use the same name as the ICON subroutines
    */

    using namespace Constants;
    namespace fm = Fast_math;
    namespace tmf = Thermo_moist_functions;

    template<typename TF> constexpr TF pi = 3.14159265359;
    template<typename TF> constexpr TF pi4 = pi<TF>/TF(4);
    template<typename TF> constexpr TF pi6 = pi<TF>/TF(6);
    template<typename TF> constexpr TF pi8 = pi<TF>/TF(8);
    template<typename TF> constexpr TF rho_w = 1.e3;                   // Density water
    template<typename TF> constexpr TF rho_0 = 1.225;                  // SB06, p48
    template<typename TF> constexpr TF pirhow = pi<TF>*rho_w<TF>/6.;
    template<typename TF> constexpr TF rho_vel = 0.4;                  // Exponent for density correction
    template<typename TF> constexpr TF kc_autocon = 9.44e+9;           // Long-Kernel
    template<typename TF> constexpr TF K_T = 2.40e-2;                  // Thermal/heat conductivity of dry air (J/m/s/K)
    template<typename TF> constexpr TF N_sc = 0.710;                   // Schmidt-Zahl (PK, S.541)
    template<typename TF> constexpr TF n_f  = 0.333;                   // Exponent von N_sc im Vent-koeff. (PK, S.541)
    template<typename TF> constexpr TF nu_l = 1.50E-5;                 // kinematic viscosity of dry air (m^2/s)

    // Limiters on ql/qr/etc.
    template<typename TF> constexpr TF q_crit = 1.e-9;                 // Min cloud/rain liquid water


    template<typename TF>
    inline TF particle_meanmass(
            Particle<TF>& particle,
            const TF q, const TF n)
    {
        // Mean mass of particle, with limiters (SB06, Eq 94)
        return std::min( std::max( q/(n+TF(Constants::dsmall)), particle.x_min ), particle.x_max );
    }

    template<typename TF>
    inline TF diffusivity(const TF T, const TF p)
    {
        // Molecular diffusivity of water vapor, from ICON.
        return TF(8.7602e-5) * exp(TF(1.81)*log(T)) / p;
    }

    template<typename TF>
    inline TF particle_diameter(
            Particle<TF>& particle,
            const TF mean_mass)
    {
        // Mass-diameter relation (SB06, Eq 32)
        return particle.a_geo * std::exp(particle.b_geo * std::log(mean_mass));
    }

    template<typename TF>
    inline TF particle_velocity(
            Particle<TF>& particle,
            const TF mean_mass)
    {
        // Mass-diameter relation (SB06, Eq 32)
        return particle.a_vel * exp(particle.b_vel * log(mean_mass)); // v = a_vel * x**b_vel
    }

    template<typename TF>
    inline TF rain_mue_dm_relation(
            Particle_rain_coeffs<TF>& coeffs,
            const TF d_m)
    {
        // mue-Dm relation of raindrops.
        TF mue;
        const TF delta = coeffs.cmu2 * (d_m - coeffs.cmu3);

        if (d_m <= coeffs.cmu3)
           mue = coeffs.cmu0 * std::tanh(fm::pow2(TF(4) * delta)) + coeffs.cmu4;
        else
           mue = coeffs.cmu1 * std::tanh(fm::pow2(delta)) + coeffs.cmu4;

        return mue;
    }

    template<typename TF>
    void setup_cloud_autoconversion(
            Particle<TF>& cloud,
            Particle_cloud_coeffs<TF>& cloud_coeffs)
    {
        const TF nu = cloud.nu;
        const TF mu = cloud.mu;

        if (mu == 1.0)  // NOTE BvS: is this safe?
        {
            //.. see SB2001
            cloud_coeffs.k_au =
                    kc_autocon<TF> / cloud.x_max * (TF(1) / TF(20)) * (nu + TF(2)) * (nu + TF(4)) / fm::pow2(nu + TF(1));
            cloud_coeffs.k_sc = kc_autocon<TF> * (nu + TF(2)) / (nu + TF(1));
        }
        else
        {
            throw std::runtime_error("Non SB01 autoconversion/selfcollection constants not (yet) setup...");
            //!.. see Eq. (3.44) of Seifert (2002)
            //cloud_coeffs%k_au = kc_autocon / cloud%x_max * (1.0_wp / 20.0_wp)  &
            //     & * ( 2.0_wp * gfct((nu+4.0)/mu)**1                           &
            //     &            * gfct((nu+2.0)/mu)**1 * gfct((nu+1.0)/mu)**2    &
            //     &   - 1.0_wp * gfct((nu+3.0)/mu)**2 * gfct((nu+1.0)/mu)**2 )  &
            //     &   / gfct((nu+2.0)/mu)**4
            //cloud_coeffs%k_sc = kc_autocon * cloud_coeffs%c_z
        }
    }

    template<typename TF>
    void diagnose_qv(
            TF* const restrict qv,
            const TF* const restrict qt,
            const TF* const restrict ql,
            const TF* const restrict qi,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;
                const int ijk = i + j*jstride + k*kstride;

                qv[ij] = qt[ijk] - ql[ij] - qi[ij];
            }
    }

    template<typename TF>
    void sedi_vel_rain(
            TF* const restrict vq,
            TF* const restrict vn,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict ql,
            const TF* const restrict rho,
            Particle<TF>& rain,
            Particle_rain_coeffs<TF>& coeffs,
            const TF rho_corr,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k,
            bool qc_present)
    {
        TF mue;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j * jstride;
                const int ijk = i + j * jstride + k * kstride;

                if (qr[ij] > q_crit<TF>)
                {
                    const TF x = particle_meanmass(rain, qr[ij], nr[ij]);
                    const TF d_m = particle_diameter(rain, x);

                    if (qc_present)
                    {
                        if (ql[ijk] >= q_crit<TF>)
                            mue = (rain.nu + TF(1.0)) / rain.b_geo - TF(1.0);
                        else
                            mue = rain_mue_dm_relation(coeffs, d_m);
                    }
                    else
                        mue = rain_mue_dm_relation(coeffs, d_m);

                    const TF d_p =
                            d_m * std::exp(TF(-1. / 3.) * std::log((mue + TF(3)) * (mue + TF(2)) * (mue + TF(1))));
                    vn[ij] = coeffs.alfa - coeffs.beta * std::exp(-(mue + TF(1)) * std::log(TF(1) + coeffs.gama * d_p));
                    vq[ij] = coeffs.alfa - coeffs.beta * std::exp(-(mue + TF(4)) * std::log(TF(1) + coeffs.gama * d_p));

                    vn[ij] *= rho_corr;
                    vq[ij] *= rho_corr;
                }
                else
                {
                    vn[ij] = TF(0);
                    vq[ij] = TF(0);
                }
            }
    }


    template<typename TF>
    void autoconversionSB(
            TF* const restrict qrt,
            TF* const restrict nrt,
            TF* const restrict qtt_liq,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict qc,
            Particle_cloud_coeffs<TF>& cloud_coeffs,
            Particle<TF>& cloud,
            Particle<TF>& rain,
            const TF cloud_rho_v,  // cloud%rho_v(i,k) in ICON, constant per layer in uHH.
            const TF Nc0,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        /*
           Autoconversion of Seifert and Beheng (2001, Atmos. Res.)
           Formation of rain by coagulating cloud droplets
        */
        const TF k_1  = 6.00e+2;   //..Parameter for Phi
        const TF k_2  = 0.68e+0;   //..Parameter fof Phi
        const TF eps  = 1.00e-25;  // NOTE BvS: dangerous for float version MicroHH.
        const TF x_s_i = TF(1) / cloud.x_max;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j * jstride;

                if (qc[ij] > q_crit<TF>)
                {
                    const TF n_c = Nc0; // cloud%n(i,k) in ICON
                    const TF x_c = particle_meanmass(cloud, qc[ij], n_c);

                    TF au = cloud_coeffs.k_au * fm::pow2(qc[ij]) * fm::pow2(x_c) * cloud_rho_v;    // NOTE `*dt` in ICON..
                    const TF tau = std::min(std::max(TF(1) - qc[ij] / (qc[ij] + qr[ij] + eps), eps), TF(0.9));
                    const TF phi = k_1 * std::pow(tau, k_2) * fm::pow3(TF(1) - std::pow(tau, k_2));
                    au = au * (TF(1) + phi / fm::pow2(TF(1) - tau));

                    nrt[ij] += au * x_s_i;
                    qrt[ij] += au;
                    qtt_liq[ij] -= au;

                    //au  = MAX(MIN(q_c,au),0.0_wp)
                    //sc  = cloud_coeffs%k_sc * q_c**2 * dt * cloud%rho_v(i,k)
                    //rain%n(i,k)  = rain%n(i,k)  + au * x_s_i
                    //rain%q(i,k)  = rain%q(i,k)  + au
                    //cloud%n(i,k) = cloud%n(i,k) - MIN(n_c,sc)
                    //cloud%q(i,k) = cloud%q(i,k) - au
                }
            }
    }


    template<typename TF>
    void accretionSB(
            TF* const restrict qrt,
            TF* const restrict qtt_liq,
            const TF* const restrict qr,
            const TF* const restrict qc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        /*
           Accretion: growth of raindrops collecting cloud droplets
           Accretion of Seifert and Beheng (2001, Atmos. Res.)
        */
        const TF k_r = 5.78;       // kernel
        const TF k_1 = 5.00e-04;   // Phi function
        const TF eps = 1.00e-25;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j * jstride;

                if (qc[ij] > TF(0) && qr[ij] > TF(0))
                {

                    // ..accretion rate of SB2001
                    const TF tau = std::min(std::max(TF(1) - qc[ij] / (qc[ij] + qr[ij] + eps), eps), TF(1));
                    const TF phi = fm::pow4(tau/(tau+k_1));
                    const TF ac  = k_r *  qc[ij] * qr[ij] * phi;  // NOTE: `*dt` in ICON..

                    qrt[ij] += ac;
                    qtt_liq[ij] -= ac;

                    //ac = MIN(q_c,ac)
                    //x_c = particle_meanmass(cloud, q_c,n_c)
                    //rain%q(i,k)  = rain%q(i,k)  + ac
                    //cloud%q(i,k) = cloud%q(i,k) - ac
                    //cloud%n(i,k) = cloud%n(i,k) - MIN(n_c,ac/x_c)
                }
            }
    }


    template<typename TF>
    void rain_selfcollectionSB(
            TF* const restrict nrt,
            const TF* const restrict qr,
            const TF* const restrict nr,
            Particle<TF>& rain,
            const TF rain_rho_v,  // rain%rho_v(i,k) in ICON, constant per layer in uHH.
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        /*
           Selfcollection & breakup: growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
           Selfcollection of Seifert and Beheng (2001, Atmos. Res.)
        */
        // Parameters based on Seifert (2008, JAS)
        const TF D_br = 1.10e-3;
        const TF k_rr = 4.33e+0;
        const TF k_br = 1.00e+3;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j * jstride;

                if (qr[ij] > TF(0))
                {
                    const TF xr = particle_meanmass(rain, qr[ij], nr[ij]);
                    const TF Dr = particle_diameter(rain, xr);

                    // Selfcollection as in SB2001
                    const TF sc = k_rr * nr[ij] * qr[ij] * rain_rho_v;  // `*dt` in ICON

                    // Breakup as in Seifert (2008, JAS), Eq. (A13)
                    TF br = TF(0);
                    if (Dr > TF(0.30e-3))
                        br = (k_br * (Dr - D_br) + TF(1)) * sc;

                    nrt[ij] += sc - br;
                }
            }
    }

    template<typename TF>
    void rain_evaporation(
            TF* const restrict qrt,
            TF* const restrict nrt,
            TF* const restrict qtt_liq,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict qv,
            const TF* const restrict ql,
            const TF* const restrict T,
            const TF* const restrict p,
            Particle_rain_coeffs<TF>& rain_coeffs,
            Particle<TF>& cloud,
            Particle<TF>& rain,
            T_cfg_2mom<TF> cfg_params,
            const TF rain_gfak,
            const TF rain_rho_v,  // rain%rho_v(i,k) in ICON, constant per layer in uHH.
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        /*
           Evaporation of rain based on Seifert (2008, J. Atmos. Sci.)
        */
        const TF D_br = TF(1.1e-3);
        const bool reduce_evaporation = true;

        const TF aa = rain_coeffs.alfa;
        const TF bb = rain_coeffs.beta;
        const TF cc = rain_coeffs.gama;

        const TF Lv2 = fm::pow2(Constants::Lv<TF>);

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j * jstride;

                const TF e_d = qv[ij] * Constants::Rd<TF> * T[ij];
                const TF e_sw = tmf::esat_liq(T[ij]); // in ICON, this calls `sat_pres_water`
                const TF s_sw = e_d / e_sw - TF(1);

                if (s_sw < TF(0) && qr[ij] > TF(0) && ql[ij] < q_crit<TF>)
                {
                    const TF D_vtp = diffusivity(T[ij], p[k]);

                    // Note that 2*pi is correct, because c_r = 1/2 is assumed
                    const TF g_d =
                        TF(2) * pi<TF> / ( Lv2 / (K_T<TF> * Constants::Rd<TF> * fm::pow2(T[ij]))
                            + Constants::Rd<TF> * T[ij] / (D_vtp * e_sw) );

                    const TF x_r = particle_meanmass(rain, qr[ij], nr[ij]);
                    const TF D_m = particle_diameter(rain, x_r);

                    // Eq. (20) of Seifert (2008)
                    TF mue;
                    if (D_m <= rain_coeffs.cmu3)
                       mue = rain_coeffs.cmu0 * tanh(
                        pow(TF(4) * rain_coeffs.cmu2 * (D_m - rain_coeffs.cmu3), rain_coeffs.cmu5)) + rain_coeffs.cmu4;
                    else
                        mue = rain_coeffs.cmu1 * tanh(
                        pow(TF(1) * rain_coeffs.cmu2 * (D_m - rain_coeffs.cmu3), rain_coeffs.cmu5)) + rain_coeffs.cmu4;

                    // Eq. (A8)
                    const TF lam = exp(TF(1)/TF(3) * log(pi6<TF> * rho_w<TF>*(mue+3.0)*(mue+2.0)*(mue+1.0)/x_r));

                    // Chebyshev approximation of Gamma(mue+5/2)/Gamma(mue+2)
                    const TF gfak =
                                   TF(0.1357940435E+01)
                        + mue * ( +TF(0.3033273220E+00)
                        + mue * ( -TF(0.1299313363E-01)
                        + mue * ( +TF(0.4002257774E-03)
                        - mue *    TF(0.4856703981E-05) ) ) );

                    const TF mm = mue + TF(5)/TF(2);

                    // Eq. (A7) rewritten with (A5) and (A9)
                    const TF f_v  =
                        rain.a_ven + rain.b_ven * pow(N_sc<TF>, n_f<TF>) * gfak *
                        sqrt(aa/nu_l<TF> * rain_rho_v / lam) *
                        (TF(1) - TF(1)/TF(2)   * (bb/aa)         * exp(mm*log(lam/(TF(1)*cc+lam)))
                               - TF(1)/TF(8)   * fm::pow2(bb/aa) * exp(mm*log(lam/(TF(2)*cc+lam)))
                               - TF(1)/TF(16)  * fm::pow3(bb/aa) * exp(mm*log(lam/(TF(3)*cc+lam)))
                               - TF(5)/TF(127) * fm::pow4(bb/aa) * exp(mm*log(lam/(TF(4)*cc+lam))) );

                    TF gamma_eva;
                    if (rain_gfak > TF(0))
                       gamma_eva = rain_gfak * (D_br/D_m) * exp(-TF(0.2)*mue); // Eq. (23)
                    else
                       gamma_eva = TF(1);

                    // Eq (A5) with (A9) and Gamma(mue+2) from (A7)
                    TF eva_q = g_d * nr[ij] * (mue + TF(1.0)) / lam * f_v * s_sw; // * dt in ICON

                    // UB: empirical correction factor to reduce evaporation of drizzle-like rain (D_m < D_br):
                    if (reduce_evaporation)
                    {
                        const TF eva_q_fak =
                            std::min(
                                std::max(
                                    (cfg_params.eva_q_fak_low +
                                    (cfg_params.eva_q_fak_high - cfg_params.eva_q_fak_low) /
                                    (cfg_params.eva_q_fak_Dbr_maxfak * D_br -
                                     cfg_params.eva_q_fak_Dbr_minfak * D_br) *
                                    (D_m - cfg_params.eva_q_fak_Dbr_minfak * D_br)),
                                cfg_params.eva_q_fak_low),
                            cfg_params.eva_q_fak_high);
                        eva_q *= eva_q_fak;
                    }

                    // Note to self: `eva_q` is a negative number.
                    // Sign diff compared to ICON is caused by their `eva_q = MAX(-eva_q,0.0_wp)`.
                    qrt[ij] += eva_q;
                    nrt[ij] += gamma_eva * eva_q / x_r;
                    qtt_liq[ij] -= eva_q;


                    //const TF eva_q = MAX(-eva_q,0.0_wp)
                    //const TF eva_n = MAX(-eva_n,0.0_wp)

                    //const TF eva_q = MIN(eva_q,q_r)
                    //const TF eva_n = MIN(eva_n,n_r)

                    //const TF atmo%qv(i,k)     = atmo%qv(i,k)     + eva_q
                    //const TF rain%q(i,k) = rain%q(i,k) - eva_q
                    //const TF rain%n(i,k) = rain%n(i,k) - eva_n
                }
            }
    }

    template<typename TF>
    void particle_particle_collection(
            TF* const restrict qpt,
            TF* const restrict qit,
            TF* const restrict nit,
            TF* const restrict qtt_ice,
            const TF* const restrict qi,
            const TF* const restrict qp,
            const TF* const restrict ni,
            const TF* const restrict np,
            const TF* const restrict Ta,
            Particle_frozen<TF>& itype,
            Particle_frozen<TF>& ptype,
            Collection_coeffs<TF>& coeffs,
            const TF rho_v,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        /*  Most simple particle-particle collection for ice particles, e.g.,
            - graupel+ice -> graupel
            - graupel+snow -> graupel
            - hail+ice -> hail
            - hail+snow -> hail
            - snow+ice -> snow

            mass from `i/itype` is collected by `p/ptype`, such that
            both the mass and concentration of `i/itype` decreases,
            and only the mass of `p/ptype` increases.
        */

        for (int j = jstart; j < jend; j++)
            #pragma ivdep
            for (int i = istart; i < iend; i++)
            {
                const int ij = i + j * jstride;

                if (qi[ij] > Sb_cold::q_crit<TF> && qp[ij] > Sb_cold::q_crit<TF>)
                {
                    //.. Sticking efficiency of Lin (1983)
                    const TF e_coll = std::min(std::exp(TF(0.09) * (Ta[ij] - Constants::T0<TF>)), TF(1));

                    const TF xp = particle_meanmass(ptype, qp[ij], np[ij]);
                    const TF dp = particle_diameter(ptype, xp);
                    const TF vp = particle_velocity(ptype, xp) * rho_v;

                    const TF xi = particle_meanmass(itype, qi[ij], ni[ij]);
                    const TF di = particle_diameter(itype, xi);
                    const TF vi = particle_velocity(itype, xi) * rho_v;

                    // Both these terms have a `* dt` in ICON; left out since we need the tendency.
                    const TF coll_n = pi4<TF> * np[ij] * ni[ij] * e_coll
                                      * (coeffs.delta_n_aa * fm::pow2(dp)
                                       + coeffs.delta_n_ab * dp * di
                                       + coeffs.delta_n_bb * fm::pow2(di))
                                  * sqrt(coeffs.theta_n_aa * fm::pow2(vp)
                                       - coeffs.theta_n_ab * vp * vi
                                       + coeffs.theta_n_bb * fm::pow2(vi)
                                       + fm::pow2(itype.s_vel));

                    const TF coll_q = pi4<TF> * np[ij] * qi[ij] * e_coll
                                      * (coeffs.delta_q_aa * fm::pow2(dp)
                                       + coeffs.delta_q_ab * dp * di
                                       + coeffs.delta_q_bb * fm::pow2(di))
                                  * sqrt(coeffs.theta_q_aa * fm::pow2(vp)
                                       - coeffs.theta_q_ab * vp * vi
                                       + coeffs.theta_q_bb * fm::pow2(vi)
                                       + fm::pow2(itype.s_vel));

                    // Gain by the destination.
                    qpt[ij] += coll_q;

                    // Loss by the source.
                    qit[ij] -= coll_q;
                    nit[ij] -= coll_n;

                    if (qtt_ice)
                        qtt_ice[ij] -= coll_q;

                    //coll_n = MIN(n_i, coll_n)
                    //coll_q = MIN(q_i, coll_q)

                    //ptype%q(i,k) = ptype%q(i,k) + coll_q
                    //itype%q(i,k)  = itype%q(i,k)  - coll_q
                    //itype%n(i,k)  = itype%n(i,k)  - coll_n
                }
            }
    }

    template<typename TF>
    void vapor_deposition_generic(
            TF* const restrict dep_q,
            const TF* const restrict prtcl_q,
            const TF* const restrict prtcl_n,
            const TF* const restrict g_i,
            const TF* const restrict s_si,
            Particle<TF>& prtcl,
            Particle_coeffs<TF>& coeffs,
            const TF rho_v,
            const TF dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {

        for (int j = jstart; j < jend; j++)
            #pragma ivdep
            for (int i = istart; i < iend; i++)
            {
                const int ij = i + j*jstride;

                if (prtcl_q[ij] == TF(0))
                    dep_q[ij] = TF(0);
                else
                {
                    const TF x = particle_meanmass(prtcl, prtcl_q[ij], prtcl_n[ij]);
                    const TF d = particle_diameter(prtcl, x);
                    const TF v = particle_velocity(prtcl, x) * rho_v;

                    TF f_v  = coeffs.a_f + coeffs.b_f * sqrt(d*v);
                    f_v = std::max(f_v, coeffs.a_f/prtcl.a_ven);

                    dep_q[ij] = g_i[ij] * prtcl_n[ij] * coeffs.c_i * d * f_v * s_si[ij] * dt;
                }
            }
    }

    template<typename TF>
    void vapor_dep_relaxation(
            // 2D Output tendencies:
            TF* const restrict qit,
            TF* const restrict nit,
            TF* const restrict qst,
            TF* const restrict nst,
            TF* const restrict qgt,
            TF* const restrict ngt,
            TF* const restrict qht,
            TF* const restrict nht,
            TF* const restrict qtt_liq,
            // Integrated deposition (so *dt):
            TF* const restrict dep_rate_ice,
            TF* const restrict dep_rate_snow,
            // 2D tmp fields:
            TF* const restrict s_si,
            TF* const restrict g_i,
            TF* const restrict dep_ice,
            TF* const restrict dep_snow,
            TF* const restrict dep_graupel,
            TF* const restrict dep_hail,
            // 2D input:
            const TF* const restrict qi,
            const TF* const restrict ni,
            const TF* const restrict qs,
            const TF* const restrict ns,
            const TF* const restrict qg,
            const TF* const restrict ng,
            const TF* const restrict qh,
            const TF* const restrict nh,
            const TF* const restrict T,
            const TF* const restrict qv,
            const TF p,
            const TF dt,
            Particle<TF>& ice,
            Particle<TF>& snow,
            Particle<TF>& graupel,
            Particle<TF>& hail,
            Particle_sphere<TF>& ice_coeffs,
            Particle_sphere<TF>& snow_coeffs,
            Particle_sphere<TF>& graupel_coeffs,
            Particle_sphere<TF>& hail_coeffs,
            const TF rho_v,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        const TF eps  = Constants::dtiny;   // 1e-20 in ICON

        // UB: if this new parameterization of n-reduction during sublimation
        //     really makes sense, move to a global constant or into the particle types.
        const bool reduce_sublimation = true;
        const TF dep_n_fac = TF(0.5);

        for (int j = jstart; j < jend; j++)
            #pragma ivdep
            for (int i = istart; i < iend; i++)
            {
                const int ij = i + j*jstride;

                if (T[ij] > Constants::T0<TF>)
                {
                    const TF e_d  = qv[ij] * Constants::Rd<TF> * T[ij];
                    const TF e_si = tmf::esat_ice(T[ij]);    // e_es(T_a) in ICON = sat_pres_ice
                    s_si[ij] = e_d / e_si - TF(1);           // Supersaturation over ice
                    const TF D_vtp = diffusivity(T[ij], p);  // D_v = 8.7602e-5 * T_a**(1.81) / p_a
                    g_i[ij] = TF(4) * Constants::pi<TF> / ( fm::pow2(Constants::Ls<TF>) /
                            (K_T<TF> * Constants::Rd<TF> * fm::pow2(T[ij])) +
                            Constants::Rd<TF> * T[ij] / (D_vtp * e_si) );
                }
                else
                {
                    g_i[ij]  = TF(0);
                    s_si[ij] = TF(0);
                }

                // Zero tmp slices
                dep_ice[ij] = TF(0);
                dep_snow[ij] = TF(0);
                dep_graupel[ij] = TF(0);
                dep_hail[ij] = TF(0);
            }

        Sb_cold::vapor_deposition_generic(
                dep_ice, qi, ni,
                g_i, s_si,
                ice, ice_coeffs,
                rho_v, dt,
                istart, iend,
                jstart, jend,
                jstride);

        Sb_cold::vapor_deposition_generic(
                dep_snow, qs, ns,
                g_i, s_si,
                snow, snow_coeffs,
                rho_v, dt,
                istart, iend,
                jstart, jend,
                jstride);

        Sb_cold::vapor_deposition_generic(
                dep_graupel, qg, ng,
                g_i, s_si,
                graupel, graupel_coeffs,
                rho_v, dt,
                istart, iend,
                jstart, jend,
                jstride);

        Sb_cold::vapor_deposition_generic(
                dep_hail, qh, nh,
                g_i, s_si,
                hail, hail_coeffs,
                rho_v, dt,
                istart, iend,
                jstart, jend,
                jstride);

        const TF zdt = TF(1) / dt;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                // Deposition only below T_3, evaporation of melting
                // particles at warmer T is treated elsewhere
                if(T[ij] > Constants::T0<TF>)
                {
                    // Depositional growth with relaxation time-scale approach based on:
                    // "A New Double-Moment Microphysics Parameterization for Application in Cloud and
                    // Climate Models. Part 1: Description" by H. Morrison, J.A.Curry, V.I. Khvorostyanov

                    const TF qvsidiff = qv[ij] - tmf::esat_ice(T[ij]) / (Constants::Rd<TF> * T[ij]);

                    if (std::abs(qvsidiff) > eps)
                    {
                        // Deposition rates are already multiplied with dt_local, therefore divide them here
                        const TF tau_i_i  = zdt / qvsidiff * dep_ice[ij];
                        const TF tau_s_i  = zdt / qvsidiff * dep_snow[ij];
                        const TF tau_g_i  = zdt / qvsidiff * dep_graupel[ij];
                        const TF tau_h_i  = zdt / qvsidiff * dep_hail[ij];

                        const TF Xi_i = ( tau_i_i + tau_s_i + tau_g_i + tau_h_i );

                        TF Xfac;
                        if (Xi_i < eps)
                            Xfac = TF(0);
                        else
                            Xfac = qvsidiff / Xi_i * (TF(1) - std::exp(-dt * Xi_i));

                        dep_ice[ij] = Xfac * tau_i_i;
                        dep_snow[ij] = Xfac * tau_s_i;
                        dep_graupel[ij] = Xfac * tau_g_i;
                        dep_hail[ij] = Xfac * tau_h_i;

                        // This limiter should not be necessary
                        if (qvsidiff < TF(0))
                        {
                            dep_ice[ij] = std::max(dep_ice[ij], -qi[ij]);
                            dep_snow[ij] = std::max(dep_snow[ij], -qs[ij]);
                            dep_graupel[ij] = std::max(dep_graupel[ij], -qg[ij]);
                            dep_hail[ij] = std::max(dep_hail[ij], -qh[ij]);
                        }

                        const TF dep_sum = dep_ice[ij] + dep_graupel[ij] + dep_snow[ij] + dep_hail[ij];

                        // NOTE BvS: this was commented out in ICON.
                        //! this limiter should not be necessary
                        //!IF (qvsidiff > 0.0_wp .and. dep_sum > qvsidiff) then
                        //!   weight = qvsidiff / dep_sum
                        //!   dep_sum          = weight * dep_sum
                        //!   dep_ice(i,k)     = weight * dep_ice(i,k)
                        //!   dep_snow(i,k)    = weight * dep_snow(i,k)
                        //!   dep_graupel(i,k) = weight * dep_graupel(i,k)
                        //!   dep_hail(i,k)    = weight * dep_hail(i,k)
                        //!END IF

                        qit[ij] += dep_ice[ij] * zdt;
                        qst[ij] += dep_snow[ij] * zdt;
                        qgt[ij] += dep_graupel[ij] * zdt;
                        qht[ij] += dep_hail[ij] * zdt;
                        qtt_liq[ij] -= dep_sum * zdt;

                        // If deposition rate is negative, parameterize the complete evaporation of some of the
                        // particles in a way that mean size is conserved times a tuning factor < 1:
                        if (reduce_sublimation)
                        {
                            const TF x_i = particle_meanmass(ice, qi[ij], ni[ij]);
                            const TF x_s = particle_meanmass(snow, qs[ij], ns[ij]);
                            const TF x_g = particle_meanmass(graupel, qg[ij], ng[ij]);
                            const TF x_h = particle_meanmass(hail, qh[ij], nh[ij]);

                            const TF dep_ice_n = std::min(dep_ice[ij], TF(0)) / x_i;
                            const TF dep_snow_n = std::min(dep_snow[ij], TF(0)) / x_s;
                            const TF dep_graupel_n = std::min(dep_graupel[ij], TF(0)) / x_g;
                            const TF dep_hail_n = std::min(dep_hail[ij], TF(0)) / x_h;

                            nit[ij] += dep_n_fac * dep_ice_n * zdt;
                            nst[ij] += dep_n_fac * dep_snow_n * zdt;
                            ngt[ij] += dep_n_fac * dep_graupel_n * zdt;
                            nht[ij] += dep_n_fac * dep_hail_n * zdt;
                        }

                        dep_rate_ice[ij] += dep_ice[ij];
                        dep_rate_snow[ij] += dep_snow[ij];
                    }
                }
            }
    }




}