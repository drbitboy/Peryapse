#!/usr/bin/env python
"""
Purpose:

  Find spacecraft periapsides in a SPICE SP-Kernel (SPK)

Usage:

  python peryapse.py SPK[ SPK[ ...]] [targets] [reffrms] [--allow-center-sun]

  - At least one SPK (filename) should be provided.
  - If no targets are supplied, only spacecraft targets will be checked.
  - [reffrms] are reference frames e.g. ECLIPJ2000; if none are supplied
    then only the frame of each segment (ephemeris) will be used in the
    output for that segment.
  - [--allow-center-sun]:  default is to ignore targets segments where
                           the either Sun or SSB is the center; this option
                           overrides that default.

Example:

    wget -q https://naif.jpl.nasa.gov/pub/naif/DEEPIMPACT/kernels/spk/zzarchive/spk_drm228_WithBurn-full.bsp

    python   peryapse.py   spk_drm230_WithBurn-full.bsp

  - Results:

    - One spacecraft (negative SPICE integer ID) was found in the SPK
      - Deep Impact Flyby Spacecraft; SPICE ID = -140, a.k.a. DIF
    - There are five segments (time periods with ephemerides) of DIF
      - Three of those DIF segments have the Sun as the center
        - Sun-centered segments are ignored by default
      - One segment has Hartley 2 as the center
      - One segment has Earth as the center
    - The Hartley 2- and Earth-centered segments are searched for periapsides

    - The results are below, with ***comment lines*** and blank lines
      inserted here for clarity

    ***Provenance of Hartley 2-centered segment; segment uses J2000 frame***

    {'SPK': 'spk_drm230_WithBurn-full.bsp',
     'calets': ['2010 NOV 03 14:13:10.156', '2010 NOV 05 13:34:13.564'],
     'center_name': 'HARTLEY 2',
     'ets': array([3.42065590e+08, 3.42236054e+08]),
     'spk_reffrm_name': 'J2000',
     'target_name': 'DEEP IMPACT FLYBY SPACECRAFT'}

    ***One Hartley 2-centered periapse found; summary in segment frame:***

    {'ET': 342150821.86116034,
     'ETCAL': '2010 NOV 04 13:53:41.861',
     'ETDIFF': 5.960464477539062e-07,
     'GFPOS': [332.52857055786126, -269.62266871771953, -553.946464467203],
     'GFVEL': [6.376120381656822, -7.452886630380785, 7.455069495989343],
     'PERPET': 342150821.86116093,
     'PERPOS': [332.52857435832516, -269.62267315998616, -553.9464600236353],
     'PERVEL': [6.376120381656822, -7.452886630380785, 7.455069495989343],
     'POSDIFF': 7.343197718544719e-06,
     'REFFRM': 'J2000'}

    ***Provenance of Earth-centered segment:***

    {'SPK': 'spk_drm230_WithBurn-full.bsp',
     'calets': ['2010 JUN 25 00:07:36.045', '2010 JUN 30 19:56:52.022'],
     'center_name': 'EARTH',
     'ets': array([3.30696456e+08, 3.31199812e+08]),
     'spk_reffrm_name': 'J2000',
     'target_name': 'DEEP IMPACT FLYBY SPACECRAFT'}

    ***One Earth-centered periapse found; summary in segment frame:***

    {'ET': 330948293.7000424,
     'ETCAL': '2010 JUN 27 22:04:53.700',
     'ETDIFF': 1.1920928955078125e-06,
     'GFPOS': [-13153.420113938717, -19465.07483123446, -28402.326001969726],
     'GFVEL': [0.6270029702688298, 4.556360226784571, -3.4129995635585164],
     'PERPET': 330948293.7000436,
     'PERPOS': [-13153.42011319127, -19465.074825802854, -28402.32600603834],
     'PERVEL': [0.6270029703936334, 4.5563602269692485, -3.4129995632890098],
     'POSDIFF': 6.827489834866299e-06,
     'REFFRM': 'J2000'}

Legend of summaries:

- ET:  ET of Geometry Finder solution, s past J2000 epoch
- ETCAL:  Leapsecond-free Calendar time of ET
- REFFRM:  Reference frame of positions and velocities
- GFPOS:  Position of target wrt center at ET, km
- GFVEL:  Velocity of target wrt center at ET, km/s
- PERPET:  ET of nearby w/position perpendicular to velocity
- PERPOS:  Position of target wrt center at PERPET, km
- PERVEL:  Velocity of target wrt center at PERPET, km/s
- ETDIFF:  GF vs. PERP time difference PERPET-ET, s
- ETDIFF:  GF vs. PERP position difference magnitude, km

N.B. Positions and velocities do include corrections for neither
     light-time nor stellar aberration (abcorr="NONE")

"""
import os
import sys
import pprint
import spiceypy as sp
import traceback as tb

do_debug = 'DEBUG' in os.environ


########################################################################
class TARGETS(object):
  """
Keep track of targets of interest:  either specifically requested
targets, or, if none were requested, then a spacecraft

"""

  def __init__(self):
    self.dt_ids = dict()
    self.deny_center_sun()

  def deny_center_sun(self): self.center_sun_allowed = False
  def allow_center_sun(self): self.center_sun_allowed = True

  def add(self,target_name):
    """Add SPICE body name as integer"""
    target_id = sp.bods2c(target_name)
    self.dt_ids[target_id] = target_name

  def of_interest(self,target_id,center_id):
    """Determine if target+center IDs are of interest;
return names or False,False

"""
    ### Ignore Sun (10) and Solar System Barycenter (0) unless
    ### asked not to
    if (10!=center_id and 0!=center_id) or self.center_sun_allowed:


      ### If no specific targets, then any negative target ID
      ### (spacecraft) is of interest
      if not self.dt_ids and target_id < 0:
        return sp.bodc2s(target_id),sp.bodc2s(center_id)

      ### Return target ID if it is in specific targets dict
      if target_id in self.dt_ids:
        return self.dt_ids[target_id],sp.bodc2s(center_id)

    return False,False


########################################################################
def do_main(lt_argv):


  tsoi,lt_reffrms,lt_spks = TARGETS(),list(),list()

  for arg in lt_argv:

    if '--allow-center-sun' == arg:
      tsoi.allow_center_sun()
      continue

    ### Try to load this arg as an SPK
    try:
      handle = sp.spklef(arg)
      sp.spkuef(handle)
      lt_spks.append(arg)
      continue
    except sp.utils.support_types.SpiceyError:
      if do_debug:
        tb.print_exc()
        print('\nIgnored spiceypy error\n')
    except:
      raise

    ### Check if this arg is an integer or name of a known inertial frame
    try:
      frmid = sp.namfrm(arg)
      if frmid != 0:
        lt_reffrms.append(sp.frmnam(frmid))
        continue
    except sp.utils.support_types.SpiceyError:
      if do_debug:
        tb.print_exc()
        print('\nIgnored spiceypy error\n')
    except:
      raise

    ### Check if this is a integer or name of a known spacecraft
    try:
      tsoi.add(arg)
      continue
    except sp.utils.support_types.SpiceyError:
      if do_debug:
        tb.print_exc()
        print('\nIgnored spiceypy error\n')
    except:
      raise

    ### anything else is a SPICE kernel e.g. TK with name/ID translation
    sp.furnsh(arg)

  ### Loop over SPKs
  for spkfn in lt_spks:

    ### Load the SPK and start a segment array search
    handle = sp.spklef(spkfn)
    sp.dafbfs(handle)

    lt_found = list()

    ### Loop over segments
    while sp.daffna():

      ### Get and parse segment summary
      ets,ic = sp.dafus(sp.dafgs(),2,6)
      (target_id,center_id,spk_reffrm_id,segtype,addrinit,addr_final
      ,)= map(int,ic)

      ### Get target name and center names
      ### - Method tsoi.of_interest will return False,False
      ###   for segments that should be ignored
      target_name,center_name = tsoi.of_interest(target_id,center_id)

      if target_name:
        ### Add ETs, names and frame is segment is of interest
        lt_found.append((ets,target_name,center_name,spk_reffrm_id,segtype,))

    ### Loop over segments of interest
    for (ets,target_name,center_name,spk_reffrm_id,segtype) in lt_found:

      ### Convert frame ID to name
      spk_reffrm_name = sp.frmnam(spk_reffrm_id)

      ### Print provenance of segment
      print()
      pprint.pprint(dict(ets=list(ets)
                        ,calets=list(map(sp.etcal,ets))
                        ,target_name=target_name
                        ,center_name=center_name
                        ,segment_type=segtype
                        ,spk_reffrm_name=spk_reffrm_name
                        ,SPK=spkfn
                        ))

      ### Initialize CNFINE argument to be passed to GFPOSC
      cnfine = sp.cell_double(2)
      sp.wninsd(ets[0],ets[1],cnfine)

      ### Maintain flag of whether a local minimum of range was found
      no_periapsis_found = True

      ### Loop over segment frame, plus any frames provided via argv
      for eval_reffrm in [spk_reffrm_name
                         ] + [reffrm for reffrm in lt_reffrms
                              if reffrm!=spk_reffrm_name
                              ]:

        ### Find any periapsides
        result = sp.cell_double(1000)
        sp.gfposc(target_name,eval_reffrm,"NONE",center_name
                 ,"RA/DEC", "RANGE","LOCMIN"
                 ,0.0,0.0,sp.spd()/4.0
                 ,2000,cnfine,result
                 )

        ### Loop over periapsides
        count = sp.wncard(result)
        for i in range(count):
          no_periapsis_found = False

          ### Get state time of periapse
          et0 = sp.wnfetd(result,i)[0]
          (gfstate,lighttime
          ,) = sp.spkezr(target_name,et0,eval_reffrm,"NONE",center_name)
          gfpos,gfvel = perppos,perpvel = gfstate[:3],gfstate[3:]

          ### Find nearby [velocity perpendicular to position] time
          perpet = et0
          deltat = - sp.vdot(gfpos,gfvel) / sp.vdot(gfvel,gfvel)
          niter = 0
          while niter < 10 and 1e-8 < abs(deltat):
            niter += 1
            perpet += deltat
            (perpstate,lighttime
            ,) = sp.spkezr(target_name,perpet,eval_reffrm,"NONE",center_name)
            perppos,perpvel = perpstate[:3],perpstate[3:]
            deltat = - sp.vdot(gfpos,gfvel) / sp.vdot(gfvel,gfvel)

          ### Print results
          pprint.pprint(dict(REFFRM=eval_reffrm
                            ,ET=et0
                            ,ETCAL=sp.etcal(et0)
                            ,GFPOS=list(gfpos)
                            ,GFVEL=list(gfvel)
                            ,ETDIFF=perpet-et0
                            ,PERPET=perpet
                            ,PERPOS=list(perppos)
                            ,PERVEL=list(perpvel)
                            ,POSDIFF=sp.vnorm(sp.vsub(gfpos,perppos))
                            ,))

      if no_periapsis_found:
        ### Flag cases where no local range minima were in the segment
        print(dict(Result='No periapsides found in this time range'))

    ### UNLOAD the SPK
    sp.spkuef(handle)


########################################################################
if "__main__" == __name__ and sys.argv[1:]:
  do_main(sys.argv[1:])