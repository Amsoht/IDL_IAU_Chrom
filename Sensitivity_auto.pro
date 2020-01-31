;-----------------------------------------------------------------------------------------------------------------------
;+
; NAME:
; SENSITIVITY_AUTO
;
; SouthTRAC ECD: auto integrate all flights for sensitivity studies
;
; AUTHOR: Thomas Wagenhäuser
;-
;------------------------------------------------------------------------------------------------------------------------
@def_ref_structs
@def_common
@wid_int_tools
@arr_get_finite_values
@conc_date
@FreeVar
@jultime2timestring
@plotsym
@strreplace
@valid_num
@tools_lib
@calc_tic
@msinfo_tools
@chrom_operations
@wid_main_tools
@wid_main
@wid_fragrat
@wid_integration
@wid_mmassviewer
@wid_mchromviewer
@wid_tpshskviewer
@wid_plotctrls
;------------------------------------------------------------------------------------------------------------------------
PRO SENSITIVITY_AUTO

  !Except = 0

  def_common, '5.19' ; string == version

  COMMON DATA
  COMMON COM_PLOT

  path = 'D:\'

  error_handler_IO = 0


  base_dir = 'C:\WiMi_ab_Maerz2019\Southtrac\Daten_Kampagnenfluege\'
  sensitivity_path = base_dir + 'Sensitivity_studies\'
  FILE_MKDIR, sensitivity_path

  chrom_dir = [$
;    '20190906b_ST04_Transfer_DE_Sal\ECD\20190906.chrom', $
;    '20190908_ST06_Transfer_Sal_BA\ECD\20190908.chrom', $
;    '20190909_ST07_Transfer_BA_RG\ECD\20190909.chrom', $
;    '20190911_ST08_Lokalflug_RG1\ECD\20190911.chrom', $
;    '20190912_Cal_ECD\ECD\20190912.chrom', $
;    '20190913_ST09_Lokalflug_RG2\ECD\20190913.chrom', $
;    '20190916_ST10_Lokalflug_RG3\ECD\20190916.chrom', $
;    '20190916a_Cal_ECD\20190916a.chrom', $
;    '20190918_ST11_Lokalflug_RG4\ECD\20190918.chrom', $
;    '20190920_ST12_Lokalflug_RG5\ECD\20190920.chrom', $
;    '20190924_ST13_Lokalflug_RG6\ECD\20190924.chrom', $
;    '20190926_Cal_ECD\20190926a.chrom', $
;    '20190926_ST14_Lokalflug_RG\ECD\20190926.chrom', $
;    '20190927_Cal_ECD\20190927.chrom', $
    '20190929_ST15_Lokalflug_RG\ECD\20190929.chrom', $
    '20190930_ST16_Lokalflug_RG\ECD\20190930.chrom', $
    '20191002_ST17_Lokalflug_RG\ECD\20191002.chrom', $
    '20191006_ST18_Transfer_RG_BA\ECD\20191006.chrom', $
    '20191007_ST19_Transfer_EZE_Sal\ECD\20191007.chrom', $
    '20191009_ST20_Transfer_Sal_EDMO\ECD\20191009.chrom', $
    '20191102_ST21_Transfer_EDMO_Sal\ECD\20191102.chrom', $
    '20191104_ST22_Transfer_Sal_EZE\ECD\20191104.chrom', $
    '20191106_ST23_Transfer_EZE_RG\ECD\20191106.chrom', $
    '20191109_ST24_Lokalflug_RG_2_1\ECD\20191109.chrom', $
    '20191112_ST25_Lokalflug_RG_2_2\ECD\20191112.chrom', $
    '20191115_ST26_Lokalflug_RG_2_3\ECD\20191115.chrom', $
    'Sensitivity_studies\20171007_10\ECD\20171007.chrom' $
    ]

;  chrom_dir = [$
;        '20190918_ST11_Lokalflug_RG4\ECD\20190918.chrom', $
;        '20190920_ST12_Lokalflug_RG5\ECD\20190920.chrom'];, $
;            '20190927_Cal_ECD\20190927.chrom'];, $
;        '20190924_ST13_Lokalflug_RG6\ECD\20190924.chrom']

  chrom_data_path = base_dir + chrom_dir
  flight_path = FILE_DIRNAME(chrom_data_path, /MARK_DIR)
;  flight_path = LIST(strsplit(chrom_data_path, '\', /EXTRACT))
;  for i=0, (N_ELEMENTS(flight_path)-1) DO flight_path[i] = STRJOIN(flight_path[i,0:-2],'\')
;  flight_path = flight_path.toarray() + '\'
  exp_dir = [$
    'gauss',$
    'gbl',$
    'gauSGbase',$
    'SGtop',$
    'SGbase',$
    'ipol_SG'$
    ]


  FOR i=0, (N_ELEMENTS(chrom_data_path)-1) DO BEGIN ;each flight

    int_ECD, chrom_data_path[i] ;read ECD data

;from wid_main.pro: 'integration':
    FOR j=0, (N_ELEMENTS(exp_dir)-1) DO BEGIN ;each integration method per flight
      exp_path = flight_path[i] + exp_dir[j]
      msinfo_path = sensitivity_path + exp_dir[j] + '_def_ms.info'
      ;read msinfo
      IF WHERE(STRMATCH(TAG_NAMES(chrom), 'subst', /FOLD_CASE) EQ 1) EQ -1 THEN BEGIN ; msinfo not loaded yet
        refs=read_subst(DEF_FILE=msinfo_path)
        IF STRLEN(refs[0].name) EQ 0 THEN RETURN ; no substances in msinfo
        refi=create_refi()
        subst=add_ires2subst(refs, refi)
        chrom=add_subst2chrom(chrom, subst)
      ENDIF ELSE BEGIN ;from wid_integration.pro: 'reload_msinfo'
        refs=read_subst(DEF_FILE=msinfo_path)
        refi=create_refi()
        subst=add_ires2subst(refs, refi)                          ; reload subst (overwrite)
        empty_chrom = create_empty_chromstrct(chrom, /CHROM_ONLY) ; reload msinfo (overwrite)
        STRUCT_ASSIGN, chrom, empty_chrom
        chrom=empty_chrom
        chrom=add_subst2chrom(chrom, subst)
        sel_chrom = 0
        sel_name = 0
      ENDELSE


;from wid_integration.pro: bat_unlist = ['bat_int', 'bat_int_noise']
      n_chrom = N_ELEMENTS(chrom.fname)
      n_subst = N_ELEMENTS(chrom[0].subst.name)

      IF (SIZE(subst, /TYPE) NE 8 OR STRLEN(subst[0].name) EQ 0) THEN BEGIN
        MSG=DIALOG_MESSAGE('Please load defaults first.', /ERROR)
        RETURN
      ENDIF

      ;'running batch method...'
      FOR sel_name=0, n_subst-1 DO BEGIN ;each substance
        FOR sel_chrom=0, n_chrom-1 DO BEGIN ;each chromatogramm
          chrom[sel_chrom].subst[sel_name].rt_win=subst[sel_name].rt_win
          chrom[sel_chrom].subst[sel_name].bl_type=subst[sel_name].bl_type
          chrom[sel_chrom].subst[sel_name].sigma=subst[sel_name].sigma
          chrom[sel_chrom].subst[sel_name].quant=subst[sel_name].quant
          chrom[sel_chrom].subst[sel_name].method=subst[sel_name].method

          call_noisecalc, sel_chrom, sel_name, event, NOISE_UNAME='noise_def_pres', /NO_WARN
          call_integration, sel_chrom, sel_name, PLOT=0
        ENDFOR
        ;save plot_intres
        FILE_MKDIR, flight_path[i]+'plot_intres\'
        plot_intres, chrom, SEL_SUBST_IX=sel_name, SAVEPLOT=flight_path[i]+'plot_intres\', CUSTOM_EXT='_'+exp_dir[j], FILE_EXT='.png', RELATIVE=1, HIDEPL=1

      ENDFOR

      sel_name = 0
      sel_chrom = 0

;save chrom and subst file ;from wid_main.pro

      IF STRLEN(chrom[0].fname)EQ 0 THEN BEGIN                                      ; chrom does not exist
        print, 'Please load Data first.'
        BREAK
      ENDIF
      IF WHERE(STRUPCASE(TAG_NAMES(chrom)) EQ 'SUBST') EQ -1 THEN BEGIN             ; chrom does not contain results information
        print, 'No Results found. Please load MSINFO fist.'
        BREAK
      ENDIF
      cdate = conc_date(chrom[0].jdate, SYSTIME(/JULIAN), cdate1=cdate1)
      savefname=exp_path + '\chrom_' + cdate + '.dat'
      ;DIALOG_PICKFILE(TITLE='Please select path and filename to save the experiment...', $
      ;  FILE='chrom_'+cdate, /OVERWRITE_PROMPT, DEFAULT_EXTENSION='dat', /WRITE, PATH=path)
      IF STRLEN(savefname[0])EQ 0 THEN BREAK                                       ; save aborted
      ;refr_status, message='saving experiment...'

      dir=FILE_DIRNAME(savefname, /MARK_DIR)                                        ; extract directory
      tmp=FILE_BASENAME(savefname)                                                  ; extract file name

      IF STRMATCH(STRMID(tmp, 0, 6), 'chrom_', /FOLD_CASE) THEN tmp=STRMID(tmp, 6)  ; delete chrom prefix if set by user (default)
      tmp=STRSPLIT(tmp, '.', /EXTRACT, COUNT=count)                                 ; split file name where a '.' is found
      IF count GT 0 THEN BEGIN                                                      ; more than one '.' found
        extention='.'+tmp[-1]                                                       ; extract extention
        fname=STRJOIN(tmp[0:-2])                                                    ; recreate file name without extention
      ENDIF

      FILE_MKDIR, dir
      save, chrom, filename=dir+'chrom_'+fname+extention
      save, subst, filename=dir+'subst_'+fname+extention
      export_subst2msinfo, subst, FNAME=dir+fname+'_def_ms.info'
;        refr_status, message='experiment saved.'





    ENDFOR
  ENDFOR
stop

END


;integrate an experiment:
PRO int_ECD, chrom_data_path
  COMMON DATA
  COMMON COM_PLOT
  ;from wid_main.pro
  BEGIN
    IF SIZE(chrom, /TYPE) NE 8 THEN chrom=create_refd()
    IF STRLEN(chrom[0].fname) NE 0 THEN BEGIN
      ;quest=DIALOG_MESSAGE('Loaded data found. Replace?', /QUESTION, /DEFAULT_NO, DIALOG_PARENT=widid.mainwid)
      ;IF quest EQ 'Yes' THEN BEGIN
        FreeVar, chrom
        FreeVar, subst
      ;ENDIF ELSE RETURN
    ENDIF
    ;destroy_wids

    ;refr_status, message='loading files...'
    chrom=read_ecd_txt(PATH=path, DEF_FILE=chrom_data_path ,VERSION=version, /LOUD)
    IF (STRLEN(chrom[0].fname) EQ 0) THEN BEGIN
      ;refr_status, message='idle'
      RETURN
    ENDIF
    ;show_list, event, chrom

    ;data_select = refresh_dsel_msel(event, chrom)
    tot_uniqm = get_uniq_mass(chrom)
    chrom.instr_type = 4 ; 0:not defined, 1:QPMS or SFMS, 2:ALMSCO_TOFMS, 3:TW_TOFMS, 4:GhostECD, 5:AED, 6:GHGGC_FID or _ECD
    ;refr_status, message='files loaded.'
  END

END
