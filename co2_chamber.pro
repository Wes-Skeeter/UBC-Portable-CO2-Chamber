;+
; name:
;   leapyear.pro
;
; purpose:
;   checks if a year is a leap-year.
;
; category:
;   date handling
;
; calling sequence:
;   result=leapyear(year)
;
; inputs:
;   year  : integer or integer array. full-digit year, e.g. 1999 or 2007
;
; output:
;   1 if leap year
;   0 if no leap year
;
; example:
;   print, leapyear(2008)
;   print, leapyear([1800,2000,2007])
;
; revision history:
;   22-dec-99 ac
;   04-may-07 ac enhanced documentation
;-

function leapyear, year

  ny=n_elements(year)
  s=intarr(ny)
  for i=0l, ny-1 do begin
    s[i]=0
    if year[i] mod 4 eq 0 then s[i]=1
    if year[i] mod 100 eq 0 and (fix(year[i]/100) mod 4) ne 0 then s[i]=0
    ;100-year rule
  endfor
  if ny eq 1 then return, s[0] else return, s

end

;+
; name:
;   months.pro
;
; purpose:
;   returns an array[13] with every first day of year (doy) for all 12 months.
;
; category:
;   date handling
;
; calling sequence:
;   result=months(year)
;
; inputs:
;   year   : year (4-digit number)
;
; output:
;   array[13], which denotes the first day of year for all months.
;   [0] = january 1st,
;   [1] = february 1st, ...
;   [11] = december 1st,
;   [12] = january of the next year.
;
; subroutines:
;   #leapyear#
;
; example
;   print, months(2002)
;
; revision history:
;   12-feb-02 ac
;   04-may-07 ac improved documentation
;-

function months, year

  if leapyear(year[0]) eq 0 then $
    months = [1,32,60,91,121,152,182,213,244,274,305,335,366] else $
    months = [1,32,61,92,122,153,183,214,245,275,306,336,367] ;leap year
  return, months

end
;+
; name:
;   doy2dat.pro
;
; purpose:
;   converts day of year (doy) and year stored in
;   logger files to a date with day and month
;
; category:
;   date handling
;
; calling sequence:
;   result=doy2dat(doy,jjjj,format=3)
;
; inputs:
;   doy    : day of year (1...366)
;   yyyy   : 4-digit year e.g. "1999" "2001"
;
; output:
;   array [month,day] or string
;
; subroutines:
;   #months#
;
; example
;   print, doy2dat(201,1999)
;   idl > [7,20] ; i.e. july-20
;   print, doy2dat(201,2000)
;   idl > [7,19] ; i.e. july-19
;   print, doy2dat([201,300],[2000,1999])
;
; revision history:
;   13-dec-01 ac
;   18-may-07 ac, enhanced documentation
;-

function doy2dat, doy, yyyy

  ni=min([n_elements(doy),n_elements(yyyy)])
  ret=intarr(2,ni)
  for i=0l, ni-1 do begin
    first_of_month=months(yyyy[i])  ;-1
    lower_mon=where(first_of_month-1 lt doy[i])
    ret[0,i]=n_elements(lower_mon)
    ret[1,i]=doy[i]-(first_of_month[ret[0,i]-1]-1)
  endfor
  return, reform(ret)

end

;+
; name:
;       met_vap_pressure
;
; purpose:
;       this function computes the vapour pressure (hpa) from air
;       temperature (degc) and relative humidity (%). the saturation
;       vapour pressure (hpa) is calculated using the magnus formula
;       (dwd p. 36, f.4.2). set u to 100% or do not specify it to
;       obtain the saturation vapour pressure for t.
;
; category:
;       meteorology
;
; calling sequence:
;       e = met_vap_pressure(t,u)
;
; input:
;       t: air temperature (degc)
;
; optional input:
;       u: relative humidity (%), default: 100%
;
; output:
;       e: vapour pressure (hpa)
;
; examples:
;       calculate the vapour pressure for u = 70% and t = 20 degc:
;
;       idl> print, #met_vap_pressure#(20,70)
;          16.393985
;
;       calculate the saturation vapour pressure for 30 degc:
;
;       idl> print, #met_vap_pressure#(30)
;          42.490755
;
; modification history:
;       written by: dieter scherer 2000 GPL
;       modified:   13-feb-2000 dis
;-

function array_processing, a0,a1,a2,a3,a4,rep_a0=rep_a0,rep_a1=rep_a1, $
  rep_a2=rep_a2,rep_a3=rep_a3,rep_a4=rep_a4

  forward_function array_processing

  ; check arguments

  n_arg = n_params()

  if n_arg lt 2 then return, 0b

  n_a0 = n_elements(a0)

  if n_a0 eq 0 then return, 0b

  ; check a1

  n_a1 = n_elements(a1)

  if n_a0 eq n_a1 then begin
    rep_a0 = a0
    rep_a1 = a1
  endif else if n_a0 eq 1 and n_a1 gt 1 then begin
    rep_a0 = replicate(a0,n_a1)
    rep_a1 = a1
  endif else if n_a0 gt 1 and n_a1 eq 1 then begin
    rep_a0 = a0
    rep_a1 = replicate(a1,n_a0)
  endif else begin
    return, 0b
  endelse

  ; check further arguments

  if n_arg ge 3 then begin
    ok = array_processing(rep_a0,a2,rep_a0=rep_a0,rep_a1=rep_a2)
    if not ok then return, 0b
    ok = array_processing(rep_a0,a1,rep_a0=rep_a0,rep_a1=rep_a1)
    if not ok then return, 0b
  endif

  if n_arg ge 4 then begin
    ok = array_processing(rep_a0,a3,rep_a0=rep_a0,rep_a1=rep_a3)
    if not ok then return, 0b
    ok = array_processing(rep_a0,a2,rep_a0=rep_a0,rep_a1=rep_a2)
    if not ok then return, 0b
    ok = array_processing(rep_a0,a1,rep_a0=rep_a0,rep_a1=rep_a1)
    if not ok then return, 0b
  endif

  if n_arg ge 5 then begin
    ok = array_processing(rep_a0,a4,rep_a0=rep_a0,rep_a1=rep_a4)
    if not ok then return, 0b
    ok = array_processing(rep_a0,a3,rep_a0=rep_a0,rep_a1=rep_a3)
    if not ok then return, 0b
    ok = array_processing(rep_a0,a2,rep_a0=rep_a0,rep_a1=rep_a2)
    if not ok then return, 0b
    ok = array_processing(rep_a0,a1,rep_a0=rep_a0,rep_a1=rep_a1)
    if not ok then return, 0b
  endif

  return, 1b

end

function vap_pressure,t,u

  ;define mask values
  mask=-9999
  maskarray=fltarr(n_elements(t))-9999

  ; saturation water vapour pressure

  t_min = -50.9   ; degc

  n_t = n_elements(t)
  n_u = n_elements(u)

  e = t - t + mask

  p = where(t lt t_min)
  if p[0] ne -1 then e[p] = mask

  c1   = 6.1078d   ; hpa
  c2_p = 17.08085d ;
  c2_n = 17.84362d ;
  c3_p = 234.175d  ; degc
  c3_n = 245.425d  ; degc

  p = where(t ge 0 and t ne mask)
  if p[0] ne -1 then e[p] = c1*exp(c2_p*t[p]/(c3_p+t[p]))

  p = where(t lt 0 and t ne mask)
  if p[0] ne -1 then e[p] = c1*exp(c2_n*t[p]/(c3_n+t[p]))

  ; u to be considered ?

  if n_params() eq 2 then begin

    ok = array_processing(e,u,rep_a0=tmp_e,rep_a1=tmp_u)
    if not ok then return, maskarray

    p = where(tmp_u ge 0 and tmp_u le 100 and tmp_e ne mask,cnt)
    if cnt gt 0 then tmp_e[p] = tmp_u[p]/100d * tmp_e[p]

    e = tmp_e

  endif

  if n_elements(e) eq 1 then e = e[0]

  return, e

end

; main event code for chamber analysis

PRO co2_chamber_event, event

  common datablock, julian, co2, soil_temp, rh, quantum
  common subset, t0, t0_delayed, t1
  common identifyers, txt, regression, iteration, iteration_redraw, done_rescaled, done_saved
  common regression_data, a0a1, RMSE, flux, molar_density, volume
  common file_handling_all, file, outfile, lun, directory
  
  chamber_type_seletor = ['"OPAQUE"','"TRANSPARENT"']
  
  CASE TAG_NAMES(event, /STRUCTURE_NAME) OF
    'WIDGET_BUTTON': BEGIN
      WIDGET_CONTROL, event.id, GET_UVALUE = event_UV
      
      ; Retrieve the Widget Window
      wDraw = WIDGET_INFO(event.top, FIND_BY_UNAME = 'DRAW')
      WIDGET_CONTROL, wDraw, GET_VALUE = graphicWin
      
      ; Retrieve the plot with the NAME
      ; provided on plot creation
      
      p = graphicWin['CO2_CHAMBER']
      p2 = graphicWin['CO2_CHAMBER_2']
      p3 = graphicWin['CO2_CHAMBER_3']
      
      base = widget_info(event.top, FIND_BY_UNAME='ROW')
      s_sample_id = widget_info(base,find_by_uname='SAMPLE_ID')
      widget_control,s_sample_id,get_value=sample_id
      s_month_id = widget_info(base,find_by_uname='SMONTH')
      widget_control,s_month_id,get_value=s_month
      s_day_id = widget_info(base,find_by_uname='SDAY')
      widget_control,s_day_id,get_value=s_day
      s_year_id = widget_info(base,find_by_uname='SYEAR')
      widget_control,s_year_id,get_value=s_year
      s_hour_id = widget_info(base,find_by_uname='SHOUR')
      widget_control,s_hour_id,get_value=s_hour
      s_minute_id = widget_info(base,find_by_uname='SMINUTE')
      widget_control,s_minute_id,get_value=s_minute
      s_second_id = widget_info(base,find_by_uname='SSECOND')
      widget_control,s_second_id,get_value=s_second
      e_hour_id = widget_info(base,find_by_uname='EHOUR')
      widget_control,e_hour_id,get_value=e_hour
      e_minute_id = widget_info(base,find_by_uname='EMINUTE')
      widget_control,e_minute_id,get_value=e_minute
      e_second_id = widget_info(base,find_by_uname='ESECOND')
      widget_control,e_second_id,get_value=e_second
      delay_id = widget_info(base,find_by_uname='DELAY')
      widget_control,delay_id,get_value=delay
      depth_id = widget_info(base,find_by_uname='DEPTH')
      widget_control,depth_id,get_value=depth
      ctype_id = widget_info(base,find_by_uname='CTYPE')
      ctype_sel = widget_info(ctype_id,/DROPLIST_SELECT)
      
      CASE event_UV OF
        'DONE': begin
        
           yesno = dialog_message('Are you sure to finish extracting samples and close the file '+file_basename(outfile)+'? You will not be able to add more samples to '+file_basename(outfile)+' later.',/question)
           if strupcase(yesno) eq 'YES' then begin
        
            close, lun
            free_lun, lun
            WIDGET_CONTROL, event.top, /DESTROY
            ;EXIT
            
           endif 
            
        end
        'READJUST': begin
        
          sjul = julday(long(s_month),long(s_day),long(s_year),long(s_hour),long(s_minute),long(s_second))
          sjul_delayed = sjul + (delay / 86400D)
          ejul = julday(long(s_month),long(s_day),long(s_year),long(e_hour),long(e_minute),long(e_second))
          diff = ejul-sjul
          margin = abs(diff)*0.2
          p.xrange = [sjul-margin,ejul+margin]
          
          t0 = where(julian ge sjul[0], t0cnt)
          t0_delayed = where(julian ge sjul_delayed[0], t0cntd)
          t1 = where(julian le ejul[0], t1cnt)
          t0 = t0[0]
          t0_delayed = t0_delayed[0]
          t1 = t1[t1cnt-1]
          
          p.title = sample_id
          if iteration_redraw gt 0 then begin
            p2.delete
            p3.delete
          endif
          p2=PLOT(julian[t0_delayed:t1], co2[t0_delayed:t1], NAME = 'CO2_CHAMBER_2', color='green', symbol='o', /sym_filled, overplot=wDraw, axis_style=1)
          p3=PLOT(julian[t0:t0_delayed], co2[t0:t0_delayed], NAME = 'CO2_CHAMBER_3', color='red', symbol='o', /sym_filled, overplot=wDraw, axis_style=1)
          mmdif = abs((max(co2[t0:t1]) - min(co2[t0:t1])) * 0.15)
          p.yrange = [min(co2[t0:t1])-mmdif,max(co2[t0:t1])+mmdif]
          iteration_redraw = iteration_redraw + 1
          
          done_rescaled = 1B
          
        end
        'CALCULATE': begin
        
          if done_rescaled eq 0B then begin
           msg = dialog_message('Please first click on "rescale" to select the proper subset for the chamber sample.')
          endif else begin
        
          if iteration gt 0 then begin
            regression.delete
            txt.delete
          endif
          
          base = widget_info(event.top, FIND_BY_UNAME='ROW')
          s_sample_id = widget_info(base,find_by_uname='SAMPLE_ID')
          
          seconds_since_start = (julian[t0_delayed:t1]-julian[t0_delayed])*86400D
          a0a1 = linfit(seconds_since_start, co2[t0_delayed:t1])
          regression=plot([julian[t0_delayed],julian[t1]],[seconds_since_start[0],$
            seconds_since_start[n_elements(seconds_since_start)-1]]*a0a1[1]+a0a1[0],$
            color='blue',linestyle='-', overplot=wDraw, axis_style = 1)
            
          modelled_co2 = seconds_since_start*a0a1[1]+a0a1[0]
          actual_co2 = co2[t0_delayed:t1]
          soil_temperature = soil_temp[t0_delayed:t1]
          RMSE = mean(sqrt((modelled_co2 - actual_co2)^2))
          
          e = vap_pressure(15,mean(rh))
          ubc_met_density, e, 1018, 273.15+15.0, rho, rho_a, rho_v
          molar_mass = 0.0289644 ; kg/mol
          molar_density = rho / molar_mass
          
          area = 0.0079 ; area in m2
          extra_volume = float(depth) * 0.01 * area
          volume = 0.0014 + extra_volume ; volume in m3
          
          flux = molar_density * volume * a0a1[1] / area
          ystrpos = 0.69
          mean_q_disp = mean(quantum[t0_delayed:t1])
          if mean_q_disp lt 0 then mean_q_disp = 0.0
          quantum_str = '!cTransparent chamber head!cPAR = '+strcompress(string(mean_q_disp,format='(f6.1)'))+' µmol $m^{-2}$ $s^{-1}$'
          if ctype_sel eq 0 then begin
           quantum_str = ''
           ystrpos = 0.76
          endif 
      
          txt = text(0.15,ystrpos,'$Efflux$ = '+strcompress(string(flux,format='(f12.2)'))+' µmol $m^{-2}$ $s^{-1}$!cRMSE = '+$
            strcompress(string(RMSE,format='(f12.1)'))+' ppm!c'+$
            'RH(max) = '+strcompress(string(max(rh[t0_delayed:t1]),format='(f4.1)'))+' %!c'+$
            'Molar density = '+strcompress(string(molar_density,format='(f4.1)'),/remove_all)+' mol $m^{-3}$'+quantum_str,/normal, color='blue')
            
          iteration = iteration + 1
          done_saved = 1B
          
          endelse
          
        end
        'WRITE': begin
        
          if done_saved eq 0B then begin
           msg = dialog_message('Please calculate the flux first before saving data.')
          endif else begin
        
          caldat, julian[t0], o_month, o_day, o_year, o_hour, o_minute, o_second
          s_datestring = '"'+string(o_year,format='(i4.4)')+'/'+string(o_month,format='(i2.2)')+'/'+string(o_day,format='(i2.2)')+'", "' $
            +string(o_hour,format='(i2.2)')+':'+string(o_minute,format='(i2.2)')+':'+string(o_second,format='(i2.2)')+'"'
            
          caldat, julian[t1], o_month, o_day, o_year, o_hour, o_minute, o_second
          e_datestring = '"'+string(o_year,format='(i4.4)')+'/'+string(o_month,format='(i2.2)')+'/'+string(o_day,format='(i2.2)')+'", "' $
            +string(o_hour,format='(i2.2)')+':'+string(o_minute,format='(i2.2)')+':'+string(o_second,format='(i2.2)')+'"'
            
          dgmsg = strupcase(dialog_message('Write current selection ('+strcompress(string(flux,format='(f12.3)'),/remove_all)+$
            ' micromol m-2 s-1) for site "'+strupcase(strcompress(sample_id,/remove_all))+'" to file '+file_basename(outfile)+'?',/question,title='Save Data...'))
            
          if dgmsg eq 'YES' then begin
          
            printf, lun, strcompress(strjoin(['"'+STRUPCASE(strcompress(sample_id,/remove_all))+'"', s_datestring, e_datestring, string(delay,format='(i4)'), $
              string(flux,format='(f12.3)'), strcompress(string(RMSE,format='(f12.1)')), $
              strcompress(string(min(co2[t0_delayed:t1]),format='(f12.2)')),$
              strcompress(string(max(co2[t0_delayed:t1]),format='(f12.2)')),$
              strcompress(string(min(rh[t0_delayed:t1]),format='(f12.2)')),$
              strcompress(string(max(rh[t0_delayed:t1]),format='(f12.2)')),$
              strcompress(string(mean(soil_temp[t0_delayed:t1]),format='(f12.2)')),$
              strcompress(string(min(soil_temp[t0_delayed:t1]),format='(f12.2)')),$
              strcompress(string(max(soil_temp[t0_delayed:t1]),format='(f12.2)')),$
              strcompress(string(volume,format='(f12.6)')),$
              strcompress(string(molar_density,format='(f12.2)')),$
              strcompress(string(mean(quantum[t0_delayed:t1]),format='(f12.2)')),$
              strcompress(string(min(quantum[t0_delayed:t1]),format='(f12.2)')),$
              strcompress(string(max(quantum[t0_delayed:t1]),format='(f12.2)')),$
              chamber_type_seletor[ctype_sel] $
              ],', '))
              
             graph_file = directory+strupcase(strcompress(sample_id,/remove_all))+'-'+string(o_year,format='(i4.4)')+string(o_month,format='(i2.2)')+string(o_day,format='(i2.2)')+string(o_hour,format='(i2.2)')+string(o_minute,format='(i2.2)')+'.png'
             p.save, graph_file  
              
          endif
          
          endelse
          
        end
        ELSE: ; do nothing
      ENDCASE
    END
    
    'WIDGET_BASE': begin
      ; Handle base resize events. Retrieve our cached padding,
      ; and our new size.
      WIDGET_CONTROL, event.id, GET_UVALUE=pad, TLB_GET_SIZE=newSize
      wDraw = WIDGET_INFO(event.top, FIND_BY_UNAME='DRAW')
      ; Change the draw widget to match the new size, minus padding.
      xy = newSize - pad
      WIDGET_CONTROL, wDraw, $
        DRAW_XSIZE=xy[0], DRAW_YSIZE=xy[1], $
        SCR_XSIZE=xy[0], SCR_YSIZE=xy[1]
    end
    
    ELSE: ; do nothing
  ENDCASE
END

pro co2_chamber

  common datablock, julian, co2, soil_temp, rh, quantum
  common identifyers, txt, regression, iteration, iteration_redraw, done_rescaled, done_saved
  common file_handling_all, file, outfile, lun, directory

  iteration = 0L
  iteration_redraw = 0L
  done_rescaled = 0B
  done_saved = 0B
  
  info = routine_info('co2_chamber',/source)
  path = file_dirname(info.path)
  ascii_def_file = path+path_sep()+'co2_chamber_ascii_def.sav'
  exists=file_search(ascii_def_file)
  if exists[0] ne '' then restore, ascii_def_file
  
  file = dialog_pickfile(path=path+path_sep(),title='Choose a data logger file from CO2 chamber system...', filter='*.dat')
  
  directory = strmid(file,0,strlen(file)-4)+'-output'+path_sep() 
  file_mkdir, directory
  outfile = directory + file_basename(file,'.dat')+'-output.csv'
  exists = file_search(outfile)
  
  if exists[0] ne '' then begin
  
    yesno = dialog_message('A file '+file_basename(outfile)+' with processed data already exists. Do you want to overwrite this file?',/question)
    if strupcase(yesno) eq 'NO' then goto, abort
  
 endif
  
  openw, lun, outfile, /get_lun
  printf, lun, '"SiteID","StartDate","StartTime","EndDate","EndTime","Delay_Seconds","CO2Flux_micromolm-2s-1","RMSE_ppm","CO2min_ppm","CO2min_max","RHmin_%","RHmax_%","SoilTavg_degC","SoilTmin_degC","SoilTmax_degC","ChamberVolume_m3","MolarDensity_molm-3","PARmean_micromolm-2s-1","PARmin_micromolm-2s-1","PARmax_micromolm-2s-1","ChamberType"'
  data = read_ascii(file, template=ascii_def)
  
  ;create proper time axis
  nt = n_elements(data.id)
  julian = dblarr(nt) * !values.f_nan
  date = double(doy2dat(data.doy,data.year))
  
  hours = reform(floor(double(data.time)/100))
  minutes = reform(double(data.time) mod 100)
  seconds = reform((data.record))
  
  julian = julday(double(reform(date[0,*])),reform(date[1,*]),reform(data.year),hours,minutes,seconds)
  caldat, julian, d_month, d_day, d_year, d_hour, d_minute, d_second
  
  s_day = string(d_day[0],format='(i2)')
  s_month = string(d_month[0],format='(i2)')
  s_year = string(d_year[0],format='(i4.4)')
  s_hour = string(d_hour[0],format='(i2)')
  s_minute = string(d_minute[0],format='(i2)')
  s_second = string(floor(d_second[0]),format='(i2)')
  
  e_day = string(d_day[nt-1],format='(i2)')
  e_month = string(d_month[nt-1],format='(i2)')
  e_year = string(d_year[nt-1],format='(i4.4)')
  e_hour = string(d_hour[nt-1],format='(i2)')
  e_minute = string(d_minute[nt-1],format='(i2)')
  e_second = string(floor(d_second[nt-1]),format='(i2)')
  def_delay = string(20,format='(i2)') ; delay for 20 seconds according tp R. Jassal
  def_depth = string(0,format='(i1)') ; delay for 20 seconds according tp R. Jassal
  def_ctype = 0B
  ;interactive plot
  
  base1 = WIDGET_BASE(/COLUMN, TITLE='CO2 Portable Chamber System - Data Analyzer - Version 1.1 (2014-07-31)', $
    /TLB_SIZE_EVENTS)
    
  wDraw = WIDGET_WINDOW(base1, UVALUE='draw', UNAME='DRAW', xsize=1160)
  
  ; Create the base for the button:
  base2 = WIDGET_BASE(base1, ROW=1, /ALIGN_LEFT, UVALUE='ROW', UNAME='ROW')
  base3 = WIDGET_BASE(base1, ROW=1, /ALIGN_LEFT, UVALUE='ROW', UNAME='ROW')
  
  
  
  ; Create the action buttons.

  label0 = WIDGET_LABEL(base2, VALUE='Site ID', xsize=65)
  sid = WIDGET_TEXT(base2, /editable, VALUE='Enter', UNAME = 'SAMPLE_ID', xsize=10)
  

  label1 = WIDGET_LABEL(base2, VALUE='Date', xsize=27)
  sday = WIDGET_TEXT(base2, /editable, VALUE=s_day, UNAME = 'SDAY', xsize=2)
  label1 = WIDGET_LABEL(base2, VALUE='/', xsize=5)
  smonth = WIDGET_TEXT(base2, /editable, VALUE=s_month, UNAME = 'SMONTH', xsize=2)
  label1 = WIDGET_LABEL(base2, VALUE='/', xsize=5)
  syear = WIDGET_TEXT(base2, /editable, VALUE=s_year, UNAME = 'SYEAR', xsize=4)
  label3 = WIDGET_LABEL(base2, VALUE='Start', xsize=40)
  shour = WIDGET_TEXT(base2, /editable, VALUE=s_hour, UNAME = 'SHOUR', xsize=2)
  label3 = WIDGET_LABEL(base2, VALUE=':', xsize=10)
  sminute = WIDGET_TEXT(base2, /editable, VALUE=s_minute, UNAME = 'SMINUTE', xsize=2)
  label3 = WIDGET_LABEL(base2, VALUE=':', xsize=10)
  esecond = WIDGET_TEXT(base2, /editable, VALUE=s_second, UNAME = 'SSECOND', xsize=2)
  label3 = WIDGET_LABEL(base2, VALUE='End', xsize=30)
  ehour = WIDGET_TEXT(base2, /editable, VALUE=e_hour, UNAME = 'EHOUR', xsize=2)
  label3 = WIDGET_LABEL(base2, VALUE=':', xsize=5)
  eminute = WIDGET_TEXT(base2, /editable, VALUE=e_minute, UNAME = 'EMINUTE', xsize=2)
  label3 = WIDGET_LABEL(base2, VALUE=':', xsize=5)
  esecond = WIDGET_TEXT(base2, /editable, VALUE=e_second, UNAME = 'ESECOND', xsize=2)

  

  
  label4 = WIDGET_LABEL(base2, VALUE='Delay (s)', xsize=65)
  delay = WIDGET_TEXT(base2, /editable, VALUE=def_delay, UNAME = 'DELAY', xsize=2)
  label3 = WIDGET_LABEL(base2, VALUE=' ', xsize=5)
  label4 = WIDGET_LABEL(base2, VALUE='Extra Depth (cm)', xsize=100)
  delay = WIDGET_TEXT(base2, /editable, VALUE=def_depth, UNAME = 'DEPTH', xsize=2)
  label3 = WIDGET_LABEL(base2, VALUE=' ', xsize=5)
  label5 = WIDGET_LABEL(base2, VALUE='Chamber Type', xsize=75)
  ctypem = WIDGET_DROPLIST(base2, VALUE=['Opaque','Transparent'], UNAME ='CTYPE')
  label6 = WIDGET_LABEL(base2, VALUE=' ', xsize=5)
  
  spc = WIDGET_LABEL(base3, VALUE=' ', xsize=65, ysize = 40)
  button1 = WIDGET_BUTTON(base3, VALUE='   Readjust   ', UVALUE = 'READJUST')
  button2 = WIDGET_BUTTON(base3, VALUE='   Calculate flux   ', UVALUE = 'CALCULATE')
  button3 = WIDGET_BUTTON(base3, VALUE='   Save to file   ', UVALUE = 'WRITE')
  done = WIDGET_BUTTON(base3, VALUE = '   Quit   ', UVALUE = 'DONE')
  
  ; Realize the widget (i.e., display it on screen).
  WIDGET_CONTROL, base1, /REALIZE
  
  ; Register the widget with the XMANAGER, leaving the IDL command
  ; line active.
  XMANAGER, 'CO2_CHAMBER', base1, /NO_BLOCK
  
  ; Cache the padding between the base and the draw
  WIDGET_CONTROL, base1, TLB_GET_SIZE=basesize
  
  xpad = basesize[0] - 640
  ypad = basesize[1] - 512
  
  WIDGET_CONTROL, base1, SET_UVALUE=[xpad,ypad]
  
  ; Retrieve the newly-created Window object.
  WIDGET_CONTROL, wDraw, GET_VALUE = graphicWin
  graphicWin.SELECT
  
  co2 = reform(data.co2_ppm)
  rh = reform(data.REL_HUM)
  soil_temp = reform(data.soil_temp)
  quantum = reform(data.unknown)
  
  ; Plot #1: In position #1 on the grid defined by LAYOUT
  
  position_set = [0.11,0.11,0.90,0.90]
  
  p=PLOT(julian, co2, NAME = 'CO2_CHAMBER', $
    LINESTYLE='-', TITLE = file_basename(file), $
    YTITLE = '$CO_2$ (ppm)', XTITLE= 'Minutes', $
    THICK=1, /CURRENT, xtickunits = 'minutes', axis_style=1, position=position_set)
    
  sjul = julian[0]
  ejul = julian[nt-1]
  diff = ejul-sjul
  
  margin = abs(diff)*0.2
  xrange_set = [sjul-margin,ejul+margin]
  p.xrange = xrange_set
  
  abort:
  
end