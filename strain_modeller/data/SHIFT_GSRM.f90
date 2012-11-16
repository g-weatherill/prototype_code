PROGRAM SHIFT_GSRM

    !Applies the SHIFT (Seismic Hazard Inferred From Tectonics)
    !hypotheses [Bird & Liu, 2007, Seismol. Res. Lett.]
    !to kinematic model GSRM (Global Strain Rate Map)
    !of Kreemer et al. [2003, Geophys. J. Int.],
    !to produce a long-term (or indefinite-term, or open-window)
    !global forecast of shallow seismicity.
    !
    !Converting each 2-D strain-rate tensor to seimicity involves
    !deciding which of the types of plate boundary studied by 
    !Bird & Kagan [2004, Bull. Seismol. Soc. Am.] is "most comparable."
    !This is done with new logic, based on the global map of deformation
    !regimes presented by Kreemer et al. [2002, Geodynamics Series, v. 30]
    !in their Plate 2.  (Tables 1 & 2 of Bird & Liu [2007] cannot be used,
    !as they implicitly assumed that all subduction zones and spreading ridges
    !and oceanic transforms would be represented by discrete faults.)
    !
    !An important approximation necessary to the conversion is the assumption
    !that the (largely elastic) strain-rates measured by GPS geodesy and represented
    !by GSRM are not very different in map pattern from long-term permanent 
    !strain rates. We consider the approximation adequate at the scale of this model.
    !
    !This program was written by Peter Bird, UCLA, in collaboration with 
    !Corne Kreemer and William Holt.  Version 1.1 of 2009.07.24

    IMPLICIT NONE
    !----------------------------------------------------------------------
    REAL, PARAMETER :: thrust_factor = 2.91
    !N.B. This is the single "fudge factor" used to adjust the overall
    !     global seismicity of this model to the correct value
    !    (taken as 21,653 m>=5.66, z<=70 km CMT earthquakes/century,
    !     which is based on 6,983 in the period 1977.01.01-2009.03.31).
    !     It represents the compounded effect of 3 corrections:
    !  (1)Systematic underestimation of seismicity of thrust faults
    !     when their slip-rates are treated as continuum strain-rates 
    !     with the continuum formula (7b) of Bird & Liu [2007], 
    !     which was based on an equivalent fault dip of 45 degrees.
    !     If all thrust faults in the world had the same dip angle,
    !     then the theoretical value of this correction factor would be
    !     (1 / ( 2 * sin(dip) * cos(dip) ) ) >= 1.
    !     For example, the interplate shear zone in the shallow part
    !     of a subduction zone dips about 14 degrees, at which angle
    !     a correction factor of 2.13 would apply.
    !  (2)Bird et al. [2009?, BSSA] found that subduction zones with 
    !     relative plate velocities V>=66 mm/a have coupled thickness and
    !     moment production rates which are 21% higher than the average
    !     for all subduction zones in their study.  Since the "subduction
    !     zone" deformation regime of Kreemer et al. [2002] is mainly 
    !     composed of these faster subduction zones, a correction factor
    !     of ~1.21 is appropriate for this effect.
    !  (3)Coupled thicknesses inferred by Bird & Kagan [2004; Table 5] 
    !     were only approximate, as they depended to some extent on the 
    !     range of years sampled.  The time window considered in their study
    !     (1977.01.01-2002.09.30) was a time of lower subduction-zone and
    !     continental-thrust activity relative to the longer window now
    !     available (1977.01.01-2009.03.31), in which global shallow (z<=70km)
    !     seismicity (m>=5.66) was higher by 5.1%, apparently due to
    !     subduction and other thrust activity that was higher by ~7.5%.
    !  Overall, thrust_factor could be estimated theoretically as about
    !     2.13 * 1.21 * 1.075 = 2.77.
    !  However, the value entered above was determined empirically, by
    !     setting the global shallow earthquake production to its value
    !     during the time window of 1977.01.01 to 2009.03.31.
    !----------------------------------------------------------------------
    CHARACTER*1 :: byte, deformation_regime_byte
    CHARACTER*1, DIMENSION(:, :), ALLOCATABLE :: regime_grid
    INTEGER :: c, c_int, c1, c1_int, c2, c2_int, columns, &
             & grid_columns, grid_rows, i, ios, j, line, m, n, &
             & r, r_int, r1, r1_int, r2, r2_int, rows, T5_column
    REAL, PARAMETER :: s_per_year = 31557600.0 ! 365.25 * 24 * 60* 60
    REAL, PARAMETER :: Earth_radius_m = 6371000.0 ! 6371 km, but expressed in m
    REAL :: c_real, c_intreal, c1_real, c1_intreal, c2_real, c2_intreal, &
          & center_normal_rate, COS_latitude_radians, &
          & dilatation_rate, dx, dy, &
          & exx, eyy, exy, &
          & e1h, e1h_OCB, e1h_ridge, e1h_transform, &
          & e2h, e2h_OCB, e2h_ridge, e2h_transform, &
          & err, err_OCB, err_ridge, err_transform, &
          & factor, fraction, fraction_to_E, fraction_to_S, &
          & grid_lat_min, grid_lat_max, grid_d_lat, &
          & grid_lon_min, grid_lon_max, grid_d_lon, &
          & intraplate_seismicity, &
          & lat, lat_max, lat_min, latitude, &
          & lon, lon_max, lon_min, longitude, &
          & r_real, r_intreal, r1_real, r1_intreal, r2_real, r2_intreal, &
          & radians_per_degree, radius_rate, second_invariant, seismicity, &
          & seismicity_OCB, seismicity_ridge, seismicity_transform, &
          & t_lon, threshold_magnitude, type_longitude, &
          & wr1c1, wr1c2, wr2c1, wr2c2
    DOUBLE PRECISION :: argument, desired_threshold_moment_Nm, G_function, &
                      & integral, integral_converted, &
                      & seismicity_at_GSRM_threshold
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE:: seismicity_grid
    !-------------------------------------------------------------------------------------
    !Definitions specific to deformation-zones input file 
    ! "tectonic_areas_bird.dat" of Kreemer et al. [2002],
    !and to its internal representation:
    INTEGER :: quadrilaterals
    CHARACTER*2, DIMENSION(:), ALLOCATABLE:: quadrilateral_byte
    REAL, DIMENSION (:, :, :), ALLOCATABLE :: quadrilateral_lonlat
    !Note that longitude range in this dataset is 0.0~360.0
    !-----------------------------------------------------------------------------------------
    !Definitions specific to GSRM input file "average_strain.dat" of Kreemer et al. [2003],
    !and to its internal (filled-out with zeros) representation GSRM_grid:
    INTEGER :: GSRM_rows, GSRM_columns ! to be computed from data below:
    REAL, PARAMETER:: GSRM_lat_min =  -79.75, GSRM_lat_max =  79.75, GSRM_d_lat = 0.5, &
                    & GSRM_lon_min = -179.70, GSRM_lon_max = 179.70, GSRM_d_lon = 0.6
                    ! in decimal degrees of East longitude & North latitude.
    REAL, PARAMETER:: strain_rate_units = 1.E-9 / 3.15576E7 ! "units of 10^-9/year"
    REAL, DIMENSION(:, :, :), ALLOCATABLE :: GSRM_grid
    !Note that longitude range in this dataset is from GSRM_lon_min to (GSRM_lon_min + 360.0).
    !-----------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------------

    !Constants to describe (uniform) seismicity of the Intraplate regions of GSRM:
    REAL :: GSRM_Intraplate_area = 4.3536E+14 ! m^2; from file "tectonic_areas.txt" created by
                                              ! program "tectonic_areas_2_dig_quadrilaterals.f90"
                                              ! (PB, 2009.07) run on input "tectonic_areas_bird.dat"
                                              ! supplied by Corner Kreemer, which was based on Plate 2
                                              ! of Kreemer et al. [2002].
    INTEGER :: GSRM_Intraplate_earthquakes = 189 ! from Global CMT, 1977.01.01-2009.03.31, 
                                                 ! for depths z <= 70 km and moment magnitudes
                                                 ! m >= 5.66 (scalar moments M >= 3.47E17 N m);
                                                 ! determined by running program "GSRM_subcatalogger.f90"
                                                 !(PB, 2009.07) on input files "tectonic_areas_bird.dat"
                                                 !(see above) and "CMT_70km_m5p66_1977-Mar2009.eqc".
    REAL :: GSRM_Intraplate_window_s = 32.25 * s_per_year ! see GSRM_Intraplate_earthquakes comment, above.
    REAL :: GSRM_Intraplate_threshold_Nm = 3.47E17 ! see GSRM_Intraplate_earthquakes comment, above.
    REAL :: GSRM_Intraplate_beta = 0.63  ! based on graph in spreadsheet "CMT_GSRM_I_tGr.xls" (PB, 2009.07).
    REAL :: GSRM_Intraplate_corner = 9.0 ! based on graph in spreadsheet "CMT_GSRM_I_tGr.xls" (PB, 2009.07).
    REAL :: GSRM_Intraplate_corner_Nm ! will be computed below
    !-------------------------------------------------------------------------------------------------------
    !****************************************************************************************
    !Definitions specific to parameters imported from Table 5 of Bird & Kagan [2004]:
    CHARACTER*6,  DIMENSION(10) :: T5_column_header
    REAL :: CMT_duration_s
    REAL, DIMENSION(10) :: T5_CMT_pure_events,        T5_CMT_pure_event_rate,   &
                         & T5_CMT_threshold_moment,                             &
                         & T5_beta,                                             &
                         & T5_corner_magnitude,       T5_corner_moment,         &
                         & T5_tGR_model_moment_rate,                            &
                         & T5_length_km, T5_length_m,                           &
                         & T5_mean_velocity_mmpa,     T5_mean_velocity_mps,     &
                         & T5_assumed_dip_degrees,    T5_assumed_dip_radians,   &
                         & T5_assumed_mu_GPa,         T5_assumed_mu_Pa,         &
                         & T5_lineIntegral_Nps,                                 &
                         & T5_coupledThickness_km,    T5_coupledThickness_m,    &
                         & T5_assumed_lithosphere_km, T5_assumed_lithosphere_m, &
                         & T5_coupling 

    !GPBT5
    !   Empirical coefficients from selected rows of Table 5 of Bird & Kagan [2004, BSSA]:
 
    ! T5_column =                       1        2        3        4        5        6        7        8        9       10
    !                                 CRB      CTF      CCB      OSR      OSR      OTF      OTF      OTF      OCB      SUB
    !                                                         normal    other     slow   medium     fast
    T5_column_header =         (/"CRB   ","CTF   ","CCB   ","OSRnor","OSRoth","OTFslo","OTFmed","OTFfas","OCB   ","SUB   "/)
    T5_CMT_pure_events =       (/   285.9,   198.5,   259.4,   424.3,    77.0,   398.0,   406.9,   376.6,   117.7,  2052.8/)
    T5_CMT_threshold_moment =  (/ 1.13E17,  3.5E17,  3.5E17, 1.13E17, 1.13E17,  2.0E17,  2.0E17,  2.0E17,  3.5E17,  3.5E17/)
    T5_beta =                  (/    0.65,    0.65,    0.62,    0.92,    0.82,    0.64,    0.65,    0.73,    0.53,    0.64/)
    T5_corner_magnitude =      (/    7.64,    8.01,    8.46,    5.86,    7.39,    8.14,    6.55,    6.63,    8.04,    9.58/)
    T5_tGR_model_moment_rate = (/ 1.67E12,  3.8E12, 1.06E13,  6.7E11,  1.9E11,  6.7E12,  9.4E11,  9.0E11,  4.6E12, 2.85E14/)
    T5_length_km =             (/  18126.,  19375.,  12516.,  61807.,  61807.,  27220.,  10351.,   6331.,  13236.,  38990./)
    T5_mean_velocity_mmpa =    (/   18.95,   21.54,   18.16,   46.40,    7.58,   20.68,   57.53,   97.11,   19.22,   61.48/)
    T5_assumed_dip_degrees =   (/     55.,     73.,     20.,     55.,     55.,     73.,     73.,     73.,     20.,     14./)
    T5_assumed_mu_GPa =        (/    27.7,    27.7,    27.7,    25.7,    25.7,    25.7,    25.7,    25.7,     49.,     49./)
    T5_lineIntegral_Nps =      (/   5.5E8,   4.4E8,   6.0E8,   5.0E9,   4.7E8,   5.2E8,   5.3E8,   5.5E8,   1.2E9, 1.58E10/)
    T5_coupledThickness_km =   (/     3.0,     8.6,     18.,    0.13,    0.40,     13.,     1.8,     1.6,     3.8,     18./)
    T5_assumed_lithosphere_km =(/      6.,     12.,     13.,      8.,      8.,     14.,     14.,     14.,     14.,     26./) 
    T5_coupling =              (/    0.50,    0.72,    1.00,   0.016,    0.05,    0.93,    0.13,    0.11,    0.27,    0.69/) 

    !  and, some readily-derived quantitites:
    radians_per_degree = 0.017453292
    CMT_duration_s = 25.7474 * s_per_year  ! 1977.01.01 through 2002.09.30 <== Do NOT update!  Must match Bird & Kagan [2004] calibration study.
    DO i = 1, 10
        T5_CMT_pure_event_rate(i) = T5_CMT_pure_events(i) / CMT_duration_s
        T5_corner_moment(i) = Moment(T5_corner_magnitude(i))
        T5_length_m(i) = T5_length_km(i) * 1000.0
        T5_mean_velocity_mps(i) = T5_mean_velocity_mmpa(i) * 0.001 / s_per_year
        T5_assumed_dip_radians(i) = T5_assumed_dip_degrees(i) * radians_per_degree
        T5_assumed_mu_Pa(i) = T5_assumed_mu_GPa(i) * 1.0E9
        T5_coupledThickness_m(i) = T5_coupledThickness_km(i) * 1000.0
        T5_assumed_lithosphere_m(i) = T5_assumed_lithosphere_km(i) * 1000.0
    END DO
    !****************************************************************************************
    GSRM_Intraplate_corner_Nm = Moment(GSRM_Intraplate_corner)

    !Print informative text on-screen:
    WRITE (*, "(' ')")
    WRITE (*, "(' PROGRAM SHIFT_GSRM:')")
    WRITE (*, "(' ')")
    WRITE (*, "(' Applies the SHIFT (Seismic Hazard Inferred From Tectonics)')")
    WRITE (*, "(' hypotheses [Bird & Liu, 2007, Seismol. Res. Lett.]')")
    WRITE (*, "(' to kinematic model GSRM (Global Strain Rate Map)')")
    WRITE (*, "(' of Kreemer et al. [2003, Geophys. J. Int.],')")
    WRITE (*, "(' to produce a long-term (or indefinite-term, or open-window)')")
    WRITE (*, "(' global forecast of shallow seismicity.')")
    WRITE (*, "(' ')")
    CALL Pause()
    WRITE (*, "(' ')")
    WRITE (*, "(' Converting each 2-D strain-rate tensor to seimicity involves')")
    WRITE (*, "(' deciding which of the types of plate boundary studied by')")
    WRITE (*, "(' Bird & Kagan [2004, Bull. Seismol. Soc. Am.] is ""most comparable.""')")
    WRITE (*, "(' This is done with new logic, based on the global map of deformation')")
    WRITE (*, "(' regimes presented by Kreemer et al. [2002, Geodynamics Series, v. 30]')")
    WRITE (*, "(' in their Plate 2.  (Tables 1 & 2 of Bird & Liu [2007] cannot be used,')")
    WRITE (*, "(' as they implicitly assumed that all subduction zones and spreading ridges')")
    WRITE (*, "(' and oceanic transforms would be represented by discrete faults.)')")
    WRITE (*, "(' ')")
    CALL Pause()
    WRITE (*, "(' ')")
    WRITE (*, "(' An important approximation necessary to the conversion is the assumption')")
    WRITE (*, "(' that the (largely elastic) strain-rates measured by GPS geodesy and represented')")
    WRITE (*, "(' by GSRM are not very different in map pattern from long-term permanent')")
    WRITE (*, "(' strain rates. We consider the approximation adequate at the scale of this model.')")
    WRITE (*, "(' ')")
    WRITE (*, "(' This program was written by Peter Bird, UCLA, in collaboration with')")
    WRITE (*, "(' Corne Kreemer and William Holt.  Version 1.1 of 2009.07.24.')")
    WRITE (*, "(' ')")
    CALL Pause()
    WRITE (*, "(' ')")
    !-------------------------------------------------------------------------------------

    !Read and memorize deformation-regime cells, from "tectonic_areas_bird.dat" on UNIT = 1:
    !N.B. This file was the basis for Plate 2 in Kreemer et al. [2002] and was provided
    !by Corne Kreemer.  Note that longitudes are in the range 0~360 East.

    WRITE (*, "(' Reading deformation regimes from ""tectonic_areas_bird.dat""...')")

    !First reading is just to count number of quadrilaterals:
    OPEN (UNIT = 1, FILE = "tectonic_areas_bird.dat", STATUS = "OLD", IOSTAT = ios)
    IF (ios /= 0) THEN
        WRITE (*, "(' ERROR: Input file tectonic_areas_bird.dat was not found.')")
        CALL Pause()
        STOP
    END IF
    line = 0
    quadrilaterals = 0 
        DefReg_indefinite: DO ! as long as there are more data in file
        READ (1, "(2X,A1)", IOSTAT = ios) byte
        line = line + 1
        IF (ios < 0) EXIT DefReg_indefinite
        IF (byte == 'C') THEN
        ELSE IF (byte == 'O') THEN
        ELSE IF (byte == 'R') THEN
        ELSE IF (byte == 'S') THEN
        ELSE
            WRITE (*, "(' ERROR: Bad byte ',A1' in line ',I6)") byte, line
            CALL Pause()
            STOP
        END IF
        DO j = 1, 5
            READ (1, *) longitude, latitude
            line = line + 1
        END DO
        quadrilaterals = quadrilaterals + 1
    END DO DefReg_indefinite
    CLOSE(UNIT = 1, DISP = "KEEP")

    ALLOCATE ( quadrilateral_byte(quadrilaterals) )
    ALLOCATE ( quadrilateral_lonlat(quadrilaterals, 5, 2) )

    !Second reading memorizes contents:
    OPEN (UNIT = 1, FILE = "tectonic_areas_bird.dat", STATUS = "OLD")
    line = 0
    quadrilaterals = 0
        DefReg_less_indefinite: DO ! as long as there are more data in file
        READ (1, "(2X,A1)", IOSTAT = ios) byte
        line = line + 1
        IF (ios < 0) EXIT DefReg_less_indefinite
        quadrilaterals = quadrilaterals + 1
        IF (byte == 'C') THEN
            quadrilateral_byte(quadrilaterals) = byte
        ELSE IF (byte == 'O') THEN
            quadrilateral_byte(quadrilaterals) = byte
        ELSE IF (byte == 'R') THEN
            quadrilateral_byte(quadrilaterals) = byte
        ELSE IF (byte == 'S') THEN
            quadrilateral_byte(quadrilaterals) = byte
        ELSE
            WRITE (*, "(' ERROR: Bad byte ',A1' in line ',I6)") byte, line
            CALL Pause()
            STOP
        END IF
        DO j = 1, 5
            READ (1, *) longitude, latitude
            line = line + 1
            IF (j == 1) THEN
                type_longitude = longitude
            ELSE
                IF (longitude - type_longitude >  350.0) longitude = longitude - 360.0
                IF (longitude - type_longitude >  350.0) longitude = longitude - 360.0
                IF (longitude - type_longitude < -350.0) longitude = longitude + 360.0
                IF (longitude - type_longitude < -350.0) longitude = longitude + 360.0
            END IF
            quadrilateral_lonlat(quadrilaterals, j, 1) = longitude
            quadrilateral_lonlat(quadrilaterals, j, 2) = latitude
        END DO
    END DO DefReg_less_indefinite
    CLOSE(UNIT = 1, DISP = "KEEP")
    !-------------------------------------------------------------------------------------

    !Read and memorize GSRM strain rates, from "average_strain.dat" on UNIT = 2:
    !N.B. This file was downloaded from http://gsrm.unavco.org/model/ in 2009.07;
    !     it seems to represent Model 1.2 of May 2004, according to associated page
    !     http://gsrm.unavco.org/intro/

    WRITE (*, "(' Reading strain-rates from ""average_strain.dat""...')")

    GSRM_rows = 1 + NINT((GSRM_lat_max - GSRM_lat_min) / GSRM_d_lat)
    GSRM_columns = 2 + NINT((GSRM_lon_max - GSRM_lon_min) / GSRM_d_lon) ! adding extra duplicated column on East
    ALLOCATE ( GSRM_grid(GSRM_rows, GSRM_columns, 3) )
    GSRM_grid = 0.0 ! whole array; note that many Intraplate tensors will NOT be changed.

    OPEN (UNIT = 2, FILE = "average_strain.dat", STATUS = "OLD", IOSTAT = ios)
    IF (ios /= 0) THEN
        WRITE (*, "(' ERROR: Input file average_strain.dat was not found.')")
        CALL Pause()
        STOP
    END IF

    READ (2, *) ! line of column headers
    READ (2, *) ! blank line
    GSRM_indefinite: DO ! as long as more input lines are available

        READ (2, *, IOSTAT = ios), lat, lon, exx, eyy, exy
        IF (ios < 0) EXIT GSRM_indefinite ! EOF encountered
        exx = exx * strain_rate_units
        eyy = eyy * strain_rate_units
        exy = exy * strain_rate_units
        i = 1 + NINT((GSRM_lat_max - lat) / GSRM_d_lat)
        j = 1 + NINT((lon - GSRM_lon_min) / GSRM_d_lon)
        GSRM_grid(i, j, 1) = exx
        GSRM_grid(i, j, 2) = eyy
        GSRM_grid(i, j, 3) = exy
        
    END DO GSRM_indefinite
    CLOSE (2, DISP = "KEEP") ! average_strain.dat

    !Copy first column of tensors to last (redundant) column, to provide wrap-around and simplify logic:
    DO i = 1, GSRM_rows
        GSRM_grid(i, GSRM_columns, 1:3) = GSRM_grid(i, 1, 1:3)
    END DO
    !-------------------------------------------------------------------------------------

    WRITE (*, "(' ==================== End of preparatory inputs =========================')")
    !================================= End of preparatory inputs =========================
    !================================= Begin production computation ======================
    WRITE (*, "(' ==================== Begin production computation ======================')")

    grid_lon_min =   0.0 ! degrees East
    grid_lon_max = 360.0 ! degrees East
    !Note that longitude range for output is the same as longitude range for quadrilaterals,
    !so no proxy longitude will be needed in comparisons.
    !However, longitude range of average_strain.dat is different, so a proxy will be used.
    grid_lat_max =  90.0 ! degrees North
    grid_lat_min = -90.0 ! degrees North
    grid_d_lon = 0.1 ! degrees
    grid_d_lat = 0.1 ! degrees

    grid_rows    = 1 + NINT((grid_lat_max - grid_lat_min) / grid_d_lat)
    grid_columns = 1 + NINT((grid_lon_max - grid_lon_min) / grid_d_lon)
    ALLOCATE ( seismicity_grid(grid_rows, grid_columns) )
    seismicity_grid = 0.0 ! must be replaced, as log(0) is negative infinity and cannot be plotted.

    WRITE (*, "(' Enter threshold (minimum) magnitude: '\)")
    READ  (*, *) threshold_magnitude
    desired_threshold_moment_Nm = Moment(threshold_magnitude)
    WRITE (*, "(' Threshold (minimum) scalar moment =',ES11.3,' N m')") desired_threshold_moment_Nm

    WRITE (*, "(' ')")
    WRITE (*, "(' Assigning deformation regime to all grid points...')")
    ALLOCATE ( regime_grid(grid_rows, grid_columns) )
    regime_grid = 'I' ! unless changed below
    DO m = 1, quadrilaterals
        deformation_regime_byte = quadrilateral_byte(m) ! to be applied to all appropriate entries in regime_grid:
        lon_min = quadrilateral_lonlat(m, 1, 1)
        lon_max = lon_min
        lat_min = quadrilateral_lonlat(m, 1, 2)
        lat_max = lat_min
        DO n = 2, 5
            lon_min = MIN(lon_min, quadrilateral_lonlat(m, n, 1))
            lon_max = MAX(lon_max, quadrilateral_lonlat(m, n, 1))
            lat_min = MIN(lat_min, quadrilateral_lonlat(m, n, 2))
            lat_max = MAX(lat_max, quadrilateral_lonlat(m, n, 2))
        END DO
        !left column limit (not worry about legal range):
        c1_real = 1.0 + (lon_min - grid_lon_min)/grid_d_lon
        c1_int = NINT(c1_real)
        c1_intreal = 1.0 * c1_int
        IF (c1_real == c1_intreal) THEN
            c1 = c1_int ! allowing an exact match on the low-longitude side
        ELSE IF (c1_real < c1_intreal) THEN
            c1 = c1_int
        ELSE ! (c1_real > c1_intreal)
            c1 = c1_int + 1 
        END IF
        !right column limit (not worry about legal range):
        c2_real = 1.0 + (lon_max - grid_lon_min)/grid_d_lon
        c2_int = NINT(c2_real)
        c2_intreal = 1.0 * c2_int
        IF (c2_real == c2_intreal) THEN
            c2 = c2_int - 1 ! not allowing an exact match on the high-longitude side
        ELSE IF (c2_real < c2_intreal) THEN
            c2 = c2_int - 1
        ELSE ! (c2_real > c2_intreal)
            c2 = c2_int 
        END IF
        DO j = c1, c2
            IF ((j >= 1).AND.(j <= grid_columns)) THEN ! j is in legal range
                !top row limit (numerically smaller; not worrying about legal range):
                r1_real = 1.0 + (grid_lat_max - lat_max)/grid_d_lat
                r1_int = NINT(r1_real)
                r1_intreal = 1.0 * r1_int
                IF (r1_real == r1_intreal) THEN
                    r1 = r1_int + 1 ! not allowing exact match on the high-latitude side
                ELSE IF (r1_real < r1_intreal) THEN
                    r1 = r1_int
                ELSE ! (r1_real > r1_intreal)
                    r1 = r1_int + 1
                END IF
                !bottom row limit (numerically larger; not worrying about legal range):
                r2_real = 1.0 + (grid_lat_max - lat_min)/grid_d_lat
                r2_int = NINT(r2_real)
                r2_intreal = 1.0 * r2_int
                IF (r2_real == r2_intreal) THEN
                    r2 = r2_int ! allowing exact match on the low-latitude side
                ELSE IF (r2_real < r2_intreal) THEN
                    r2 = r2_int - 1
                ELSE ! (r2_real > r2_intreal)
                    r2 = r2_int
                END IF
                DO i = r1, r2
                    IF ((i >= 1).AND.(i <= grid_rows)) THEN ! i is in legal range
                        !store regime information for use in next operations block:
                        regime_grid(i, j) = deformation_regime_byte
                    END IF ! i is in legal range
                END DO ! i = r1, r2
            END IF ! j is in legal range
        END DO ! j = c1, c2
    END DO ! m = 1, quadrilaterals
    DO i = 1, grid_rows
        regime_grid(i, grid_columns) = regime_grid(i, 1) ! copy 1st column to last column
    END DO
    !-------------------------------------------------------------------------------------

    !Precompute intraplate_seismicity (to save duplication, and to use a lower-limit):
    seismicity_at_GSRM_threshold = GSRM_Intraplate_earthquakes / (GSRM_Intraplate_area * 1.D0 * GSRM_Intraplate_window_s)
    argument = (GSRM_Intraplate_threshold_Nm - desired_threshold_moment_Nm) / GSRM_Intraplate_corner_Nm
    IF (argument < -100.0D0) THEN
        G_function = 0.0D0
        WRITE (*, "(' WARNING: argument underflow in EXP(argument): intraplate seismicity will be zero.')")
        CALL Pause()
    ELSE
        G_function = ((desired_threshold_moment_Nm / GSRM_Intraplate_threshold_Nm)**(-GSRM_Intraplate_beta)) * &
                     & EXP(argument)
    END IF
    intraplate_seismicity = G_function * seismicity_at_GSRM_threshold
    !N.B. This uniform value will be applied to all 'I' grid points, and will be used as a lower limit 
    !     for the seismicity of other grid points as well.

    !-------------------------------------------------------------------------------------

    DO i = 1, grid_rows ! (N to S)
        latitude = grid_lat_max - (i-1) * grid_d_lat
        WRITE (*, "(' Doing row with latitude = ',F7.3,'...')") latitude

        DO j = 1, grid_columns ! (W to E)
            longitude = grid_lon_min + (j-1) * grid_d_lon
            !proxy longitude which is different by +- 360 will often be needed:
            t_lon = longitude
            IF (t_lon < GSRM_lon_min) t_lon = t_lon + 360.0
            IF ((t_lon - GSRM_lon_min) > 360.0) t_lon = t_lon - 360.0

            !recall deformation regime, stored in previous operations block:
            deformation_regime_byte = regime_grid(i, j)
          
            IF (deformation_regime_byte == 'I') THEN
                !strain-rate will not be needed; just set it to zero for sake of tidiness:
                exx = 0.0; eyy = 0.0; exy = 0.0
            ELSE ! must look up the strain-rate

                r_real = 1.0 + (GSRM_lat_max - latitude)/GSRM_d_lat
                r_int = NINT(r_real)
                c_real = 1.0 + (t_lon - GSRM_lon_min)/GSRM_d_lon ! using proxy longitude!
                c_int = NINT(c_real)
                IF ((r_int >= 1).AND.(r_int <= GSRM_rows).AND.(c_int >= 1).AND.(c_int <= GSRM_columns)) THEN
                    exx = GSRM_grid(r_int, c_int, 1)
                    eyy = GSRM_grid(r_int, c_int, 2)
                    exy = GSRM_grid(r_int, c_int, 3)
                ELSE ! point lies outside of the GSRM grid; treat as Intraplate:
                    exx = 0.0; eyy = 0.0; exy = 0.0
                END IF ! inside or outside the GSRM grid

                !If strain-rate tensor is (0., 0., 0.), then set deformation_regime_byte == 'I':
                IF ((exx == 0.0).AND.(eyy == 0.0).AND.(exy == 0.0)) deformation_regime_byte = 'I'
            END IF ! Intraplate, or not

            !Characterize the strain-rate tensor just found:
            second_invariant = SQRT(exx**2 + eyy**2 + 2.0 *(exy**2)) ! (per Kreemer et al. [2002]; used in test plots)
            dilatation_rate = exx + eyy
            err = -dilatation_rate ! based on assumption of volume conservation in permanent, non-elastic strain
            !N.B. As strain-rates of GSRM may include large elastic parts, this is an APPROXIMATION when used here!
            center_normal_rate = (exx + eyy)/2.0
            radius_rate = SQRT((exx - center_normal_rate)**2 + exy**2)
            e1h = center_normal_rate - radius_rate ! most-compressive horizontal principal value
            e2h = center_normal_rate + radius_rate ! most-extensional horizontal principal value

            !---------- Use regime-specific logic (and subprogram GSRM_continuum_seismicity)
            !           to find the rate of earthquakes/m^2/s above desired_threshold_moment_Nm:

            IF (deformation_regime_byte == 'I') THEN ! Intraplate
                seismicity_grid(i, j) = intraplate_seismicity ! (precomputed above, outside these loops)

            ELSE IF (deformation_regime_byte == 'C') THEN ! Continental
                !This classification follows Table 2 of Bird & Liu [2007]:
                IF (err > 0.0) THEN ! thrust faulting present
                    IF (err > (0.364 * e2h)) THEN 
                        T5_column = 3 ! CCB
                        factor = thrust_factor
                    ELSE ! thrust faulting is minor compared to strike-slip
                        T5_column = 2 ! CTF
                        factor = 1.0
                    END IF
                ELSE IF (err == 0.0) THEN ! pure strike-slip
                    T5_column = 2 ! CTF
                    factor = 1.0
                ELSE ! (err < 0.0); normal faulting present
                    IF (err >= (0.364 * e1h)) THEN
                        T5_column = 2 ! CTF
                        factor = 1.0
                    ELSE ! (err < 0.364*e1h) 
                        T5_column = 1 ! CRB
                        factor = 1.0
                    END IF
                END IF
                CALL SHIFT_continuum_seismicity (desired_threshold_moment_Nm, & ! INTENT(IN)
                                               & e1h, e2h, err, T5_column, &    ! INTENT(IN)
                                               & seismicity) ! INTENT(OUT)
                seismicity = seismicity * factor ! either thrust_factor, or 1.0
                seismicity_grid(i, j) = MAX(seismicity, intraplate_seismicity)

            ELSE IF (deformation_regime_byte == 'O') THEN ! Oceanic (diffuse)
                T5_column = 9 ! OCB
                CALL SHIFT_continuum_seismicity (desired_threshold_moment_Nm, & ! INTENT(IN)
                                               & e1h, e2h, err, T5_column, &    ! INTENT(IN)
                                               & seismicity) ! INTENT(OUT)
                IF (err > 0.0) THEN ! thrust faulting present
                    IF (err > (0.364 * e2h)) THEN 
                        factor = thrust_factor
                    ELSE ! thrust faulting is minor compared to strike-slip
                        factor = 1.0
                    END IF
                ELSE IF (err == 0.0) THEN ! pure strike-slip
                    factor = 1.0
                ELSE ! (err < 0.0); normal faulting present
                    factor = 1.0
                END IF
                seismicity = seismicity * factor
                seismicity_grid(i, j) = MAX(seismicity, intraplate_seismicity)

            ELSE IF (deformation_regime_byte == 'R') THEN ! Ridge/transform
                !This classification is described in file
                !"methods_memo_2008_09_19.doc" under the heading "IV".
                !It separates a general R (ridge/transform) strain-rate tensor into 2 parts,
                !and computes seismicity separately for each part.
                IF ((e1h > 0.0).AND.(e2h > 0.0)) THEN ! spreading along both horizontal principal axes; no transforms:
                    T5_column = 4 ! OSR normal
                    CALL SHIFT_continuum_seismicity (desired_threshold_moment_Nm, & ! INTENT(IN)
                                                   & e1h, e2h, err, T5_column, &    ! INTENT(IN)
                                                   & seismicity) ! INTENT(OUT)
                    seismicity_grid(i, j) = MAX(seismicity, intraplate_seismicity)
                ELSE IF (e1h == 0.0) THEN ! simple uniaxial spreading; no transforms
                    T5_column = 4 ! OSR normal
                    CALL SHIFT_continuum_seismicity (desired_threshold_moment_Nm, & ! INTENT(IN)
                                                   & e1h, e2h, err, T5_column, &    ! INTENT(IN)
                                                   & seismicity) ! INTENT(OUT)
                    seismicity_grid(i, j) = MAX(seismicity, intraplate_seismicity)
                ELSE IF (((e1h * e2h) < 0.0).AND.((e1h + e2h) >= 0.0)) THEN ! mixture of spreading & transform:
                    e1h_ridge = 0.0 ! no shortening along ridge axis
                    e2h_ridge = e1h + e2h ! all dilatation is due to spreading
                    err_ridge = -(e1h_ridge + e2h_ridge) ! volume conservation
                    T5_column = 4 ! OSR normal
                    CALL SHIFT_continuum_seismicity (desired_threshold_moment_Nm, & ! INTENT(IN)
                                                   & e1h_ridge, e2h_ridge, err_ridge, T5_column, &    ! INTENT(IN)
                                                   & seismicity_ridge) ! INTENT(OUT)
                    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    e1h_transform = e1h ! all shortening is due to transform activity
                    e2h_transform = -e1h_transform ! no dilatation due to transforms
                    err_transform = 0.0
                    T5_column = 7 ! OTF medium
                    CALL SHIFT_continuum_seismicity (desired_threshold_moment_Nm, & ! INTENT(IN)
                                                   & e1h_transform, e2h_transform, err_transform, T5_column, &    ! INTENT(IN)
                                                   & seismicity_transform) ! INTENT(OUT)
                    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    seismicity_grid(i, j) = MAX((seismicity_ridge + seismicity_transform), intraplate_seismicity)
                ELSE IF (((e1h * e2h) < 0.0).AND.((e1h + e2h) < 0.0)) THEN ! mixture of thrusting & transform:
                    e1h_OCB = e1h + e2h ! all dilatation is due to thrusting
                    e2h_OCB = 0.0 ! no change in lengths of thrust fault traces
                    err_OCB = -(e1h_OCB + e2h_OCB) ! volume conservation
                    T5_column = 9 ! OCB
                    CALL SHIFT_continuum_seismicity (desired_threshold_moment_Nm, & ! INTENT(IN)
                                                   & e1h_OCB, e2h_OCB, err_OCB, T5_column, &    ! INTENT(IN)
                                                   & seismicity_OCB) ! INTENT(OUT)
                    seismicity_OCB = seismicity_OCB * thrust_factor
                    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    e2h_transform = e2h ! all extension is due to transform activity
                    e1h_transform = -e2h_transform ! no dilatation due to transform activity
                    err_transform = 0.0
                    T5_column = 7 ! OTF medium
                    CALL SHIFT_continuum_seismicity (desired_threshold_moment_Nm, & ! INTENT(IN)
                                                   & e1h_transform, e2h_transform, err_transform, T5_column, &    ! INTENT(IN)
                                                   & seismicity_transform) ! INTENT(OUT)
                    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    seismicity_grid(i, j) = MAX((seismicity_OCB + seismicity_transform), intraplate_seismicity)
                ELSE IF ((e1h < 0.0).AND.(e2h < 0.0)) THEN ! shortening in both directions; use OCB model:
                    T5_column = 9 ! OCB
                    CALL SHIFT_continuum_seismicity (desired_threshold_moment_Nm, & ! INTENT(IN)
                                                   & e1h, e2h, err, T5_column, &    ! INTENT(IN)
                                                   & seismicity) ! INTENT(OUT)
                    seismicity = seismicity * thrust_factor
                    seismicity_grid(i, j) = MAX(seismicity, intraplate_seismicity)
                ELSE ! (not expected to happen)
                    WRITE (*, "(' ERROR: Unforseen case; fix program logic.  e1h = ',ES10.3,', e2h = ',ES10.3)") e1h, e2h
                    CALL Pause()
                    STOP
                END IF

            ELSE IF (deformation_regime_byte == 'S') THEN ! Subduction
                T5_column = 10 ! SUB
                CALL SHIFT_continuum_seismicity (desired_threshold_moment_Nm, & ! INTENT(IN)
                                               & e1h, e2h, err, T5_column, &    ! INTENT(IN)
                                               & seismicity) ! INTENT(OUT)
                seismicity = seismicity * thrust_factor
                seismicity_grid(i, j) = MAX(seismicity, intraplate_seismicity)

            ELSE ! bad byte!
                WRITE (*, "(' ERROR: deformation_regime_byte = ',A,' at latitude ',F7.3,', longitude ',F8.3)") deformation_regime_byte, latitude, longitude
                CALL Pause()
                STOP
            END IF ! selection based on deformation_regime_byte

        END DO ! j = 1, grid_columns (W to E)
    END DO ! i = 1, grid_rows (N to S)

    !-----------------------------------------------------------------------------------

    !Output the .grd file:
    WRITE (*, "(' Writing .grd file and computing area-integral...')")
    OPEN (UNIT = 3, FILE = "SHIFT_GSRM.grd", STATUS = "NEW", IOSTAT = ios)
    WRITE (3, "(3F11.5,' = lon_min, d_lon, lon_max  Gridded long-term seismicity above magnitude ',F6.3,' in events/square-meter/s')") grid_lon_min, grid_d_lon, grid_lon_max, threshold_magnitude
    WRITE (3, "(3F11.5,' = lat_min, d_lat, lat_max  from program SHIFT_GSRM, by P. Bird, UCLA  ')") grid_lat_min, grid_d_lat, grid_lat_max
    WRITE (3, "(1P,10D10.2)") ((seismicity_grid(i, j), j = 1, grid_columns), i = 1, grid_rows)
    !Note: Line breaks in output file have no relation to the logical shape of this array.
    !      Programs that read .grd files must take this into account, by reading the whole
    !      array with a single READ (as it was written), and not row-by-row.
    !-----------------------------------------------------------------------------------
   
    !Integrate global shallow seismicity: 

   !Integrate seismicity over grid area:
    integral = 0.0D0 ! (just initializing)
    DO i = 1, grid_rows
        latitude = grid_lat_max - (i - 1) * grid_d_lat
        COS_latitude_radians = COS(latitude * radians_per_degree)
        IF ((i == 1).OR.(i == grid_rows)) THEN
            dy = 0.5D0 * grid_d_lat * radians_per_degree * Earth_radius_m
        ELSE
            dy = grid_d_lat * radians_per_degree * Earth_radius_m
        END IF
        DO j = 1, grid_columns
            IF ((j == 1).OR.(j == grid_columns)) THEN
                dx = 0.5D0 * grid_d_lon * COS_latitude_radians * radians_per_degree * Earth_radius_m
            ELSE
                dx = grid_d_lon * COS_latitude_radians * radians_per_degree * Earth_radius_m
            END IF
            integral = integral + seismicity_grid(i, j) * dx * dy
        END DO
    END DO
    integral_converted = integral * 100.0 * s_per_year ! now, earthquakes per century
    WRITE (*, "(/' Area integral of seismicity =' &
              & /' ',ES11.4,' earthquakes/s = ',ES11.4,' earthquakes/century.')") integral, integral_converted
    WRITE (3, "(  'Area integral of seismicity =' &
              &     /ES11.4,' earthquakes/s = ',ES11.4,' earthquakes/century.')") integral, integral_converted
    CLOSE (UNIT = 3, DISP = "KEEP") ! seismicity.grd file
    !-----------------------------------------------------------------------------------
    !Signal completion of job:
    !================================= End production computation ======================
    WRITE (*, "(' ==================== end production computation ======================')")
    WRITE (*, "(' Job completed.')")
    CALL Pause()

    !-------------------------------------------------------------------------------------

CONTAINS

    REAL FUNCTION Moment(magnitude)
    !   Returns scalar seismic moment in N m, per 
    !   equation (1) of Bird & Kagan [2004; Bull. Seism. Soc. Am.];
    !   originally from Hanks & Kanamori [1979; J. Geophys. Res., v. 84, p. 2348-2350].
        IMPLICIT NONE
        REAL, INTENT(IN) :: magnitude
        moment = 10.**((1.5 * magnitude) + 9.05)
    END FUNCTION Moment

    SUBROUTINE Pause ()
        IMPLICIT NONE
        CHARACTER*1 c1
        WRITE (*,"(' Press [Enter] to continue...'\)")
        READ (*,"(A)") c1
    END SUBROUTINE Pause

    SUBROUTINE SHIFT_continuum_seismicity (desired_threshold_moment_Nm, &
                                         & e1h, e2h, err, T5_column, &
                                         & seismicity)
        !Computes "seismicity" (# of earthquakes/m^2/s, including aftershocks)
        !above "desired_threshold_moment_Nm",
        !based on principal long-term strain rates
        !(e1h = most-compressive horizontal; 
        ! e2h = least-compressive horizontal;
        ! err = vertical)
        !and choice of most comparable plate boundary expressed by
        !selection of column "T5_column" in the T5_{whatever} tables in
        !the main program above.
        !The only reason that this code is isolated in a subprogram
        !is that it is called twice for grid points in tectonic
        !regime "R" of the Global Strain Rate Map model,
        !and redundant code would be needed if this were placed in-line.

        !The underlying assumptions and equations are those of the SHIFT
        !(Seismic Hazard Inferred From Tectonics)
        !method of Bird & Liu [2007, Seismol. Res. Lett.].

        IMPLICIT NONE
        !But note that arrays T5_{whatever} of main program are read: INTENT(IN).
        DOUBLE PRECISION, INTENT(IN) :: desired_threshold_moment_Nm
        REAL, INTENT(IN) :: e1h, e2h, err
        INTEGER, INTENT(IN) :: T5_column
        REAL, INTENT(OUT) :: seismicity

        REAL :: CMT_threshold_moment, corner_moment, e1_rate, e2_rate, e3_rate, &
              & M_persec_per_m2, M_persec_per_m3, &
              & seismicity_at_CMT_threshold, tGR_beta
        DOUBLE PRECISION :: argument, G_function

        !Convert principal strain rates in continuum to seismic moment rate per unit volume,
        !per equation (7b) of Bird & Liu [2007]):
        e1_rate = MIN(e1h, e2h, err) ! always negative
        e3_rate = MAX(e1h, e2h, err) ! always positive
        e2_rate = 0.0 - e1_rate - e3_rate      ! may have either sign
        IF (e2_rate < 0.0) THEN ! both e1_rate and e2_rate are negative; positive e3_rate is the unique axis:
            M_persec_per_m3 = T5_assumed_mu_Pa(T5_column) * 2.0 * e3_rate
        ELSE ! e2_rate >= 0.0 ; so both e2_rate and e3_rate postive; negative e1_rate is the unique axis:
            M_persec_per_m3 = T5_assumed_mu_Pa(T5_column) * 2.0 * (-e1_rate)
        END IF

        !Convert moment rate per unit volume to moment rate per unit area:
        M_persec_per_m2 = M_persec_per_m3 * T5_coupledThickness_m(T5_column)

        !Determine subsampling-point value of seismicity (earthquakes/m**2/s),
        !per equation (4) of Bird & Liu [2007]:
        seismicity_at_CMT_threshold = T5_CMT_pure_event_rate(T5_column) * 1.0D0 * &
                                    & M_persec_per_m2 / T5_tGR_model_moment_rate(T5_column)
                                   !Caution: Must avoid underflow, as terms are likely to
                                   !         be roughly: 2E-7 * 3E-6 / 4E12 = 1E-25.
                                   !         Therefore, 1.0D0 is introduced into product
                                   !         to force double-precision computation.

        !Correct seismicity to desired_threshold_moment_Nm,
        !per equation (5) of Bird & Liu [2007]:
        !using the tapered Gutenberg-Richter distribution
        !with coefficients from Table 5 above:
        tGR_beta = T5_beta(T5_column)
        corner_moment = T5_corner_moment(T5_column)
        CMT_threshold_moment = T5_CMT_threshold_moment(T5_column)
        !G = ((desired_threshold_moment_Nm / CMT_threshold_moment)**(-beta)) * &
        !  & EXP((CMT_threshold_moment - desired_threshold_moment_Nm) / corner_moment)
        argument = (CMT_threshold_moment - desired_threshold_moment_Nm) / corner_moment
        IF (argument < -100.0D0) THEN
            G_function = 0.0D0
        ELSE
            G_function = ((desired_threshold_moment_Nm / CMT_threshold_moment)**(-tGR_beta)) * &
                         & EXP(argument)
        END IF
        seismicity = G_function * seismicity_at_CMT_threshold

    END SUBROUTINE SHIFT_continuum_seismicity

END PROGRAM SHIFT_GSRM
