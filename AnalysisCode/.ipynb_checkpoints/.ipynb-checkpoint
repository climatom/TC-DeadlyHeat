{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# TC Met-Impact\n",
    "Here the impact of TCs on the surface meteorology is assessed. Prior to running this notebook, the WFDEI gridded fields must have been extracted to the points of TC landfall. \n",
    "\n",
    "The main analysis components covered here are: \n",
    "\n",
    "(1) The generation of Day-of-Year climatologies for each landfalling site\n",
    "\n",
    "(2) Extraction of the met time series for the 30-days post TC landfall at each site, expressed as anomalies\n",
    "\n",
    "Under (2) the maximum absolute HI per TC (rather than site) is also assessed, so that we can assess whether conditions became \"dangerously hot\". \n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 147,
   "metadata": {},
   "outputs": [],
   "source": [
    "# Perform imports and set parameters\n",
    "import GeneralFunctions as GF\n",
    "import TC_Utils as tc\n",
    "import calendar, pickle, numpy as np\n",
    "from netCDF4 import Dataset\n",
    "import matplotlib.pyplot as plt\n",
    "\n",
    "# MetCDF4 file -- holding WFDEI data (extracted to TC landfall points)\n",
    "metfile=\"/media/gytm3/WD12TB/TropicalCyclones/TC-DeadlyHeat/Data/TC_met_all.nc\"\n",
    "# has variables...\n",
    "vs=[\"hi\",\"Tair\",\"Qair\",\"PSurf\"]\n",
    "\n",
    "# Use a stdv of...\n",
    "stdv=15./1.96 # to control the \"bandwidth\" of the Gaussian kernel\n",
    "\n",
    "# Examine met this many days after TC landfall\n",
    "ndays=30\n",
    "\n",
    "# Examine this many days BEFORE TC landfall\n",
    "nbefore=30\n",
    "\n",
    "# Main meta-file with all TC landfalling info \n",
    "locfile=\"/media/gytm3/WD12TB/TropicalCyclones/TC-DeadlyHeat/Data/LandFall.txt\"\n",
    "# NOTE>> Format is:\n",
    "# [0] ID; [1] Year [2] jd.dayfrac; [3] lat; [4] lon; [5] grid lat\n",
    "# [6] grid lon; [7] grid row; [8] grid col; [9] dist between gridpoint and \n",
    "# TC landfall coordinate\n",
    "\n",
    "# Now read...\n",
    "meto=Dataset(metfile,\"r\")\n",
    "locs=np.loadtxt(locfile,skiprows=1)\n",
    "# get ntime and nlocs\n",
    "ntime,nlocs=meto.variables[vs[0]].shape\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## [1] Computing the seasonal cycle (and anomalies)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 65,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Finished with variable: hi\n",
      "Finished with variable: Tair\n",
      "Finished with variable: Qair\n",
      "Finished with variable: PSurf\n"
     ]
    }
   ],
   "source": [
    "# Prepare met data for iteration to compute the seasonal cycles\n",
    "time=meto.variables[\"time\"]\n",
    "yr,mon,day,hr,dt=\\\n",
    "GF.conTimes(time_str=time.units,calendar=time.calendar,\\\n",
    "                    times=time[:],safe=False)\n",
    "grid_times=np.floor(tc.dt2decday(dt)) # decimal time (WFDEI)\n",
    "gridlat=meto.variables[\"lat\"][:]\n",
    "gridlon=meto.variables[\"lon\"][:]\n",
    "\n",
    "# Now we're going to iterate over the met vars and the locations --\n",
    "# calling the kernel DoY smoother to compute the DoY climatology \n",
    "# at each site (DoYs 1-->366)\n",
    "clims={}; anoms={}\n",
    "for v in vs:\n",
    "    clims[v]=np.array([tc.kernel_smooth_doy(meto.variables[v][:,ii],grid_times,stdv)[1] for ii in range(nlocs)])\n",
    "    anoms[v]=np.zeros((ntime,nlocs))\n",
    "    # Compute anomalies by looping over unique days\n",
    "    for dd in range(1,367):\n",
    "        idx=grid_times==dd\n",
    "        anoms[v][idx,:]=meto.variables[v][idx,:]-clims[v][:,dd-1]\n",
    "    print \"Finished with variable: %s\" % v\n",
    "\n",
    "    \n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## [2] Extracting the anomalies before/after TC landfall"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 137,
   "metadata": {},
   "outputs": [],
   "source": [
    "\n",
    "# Now the post-TC met should be found for the ndays post TC-landfall. \n",
    "# The easiest way to do this is to use the year JD.xx of landfall and then take the year and jd+1:jd+ndays using \n",
    "# yr and grid_times. That logical index can be used to slice anoms and the raw meto.variables variable\n",
    "\n",
    "# Make substitutions for clearer code\n",
    "tc_yy=locs[:,1]\n",
    "tc_jd=np.floor(locs[:,2])\n",
    "nim=0\n",
    "for v in vs:\n",
    "    \n",
    "    impact[v]=np.zeros((nlocs,ndays))\n",
    "    profile[v]=np.zeros((nlocs,(ndays+nbefore+1)))\n",
    "    for ll in range(nlocs):\n",
    "        \n",
    "        \n",
    "        if tc_yy[ll]>2015: continue # we do this because only HI goes up to 2016/17 -- was extended with ERA-I. \n",
    "            \n",
    "        if tc_yy[ll]==2015 and tc_jd[ll]>(365-ndays): continue # again, because 30-day window will trigger an \n",
    "            # out-of-bounds error\n",
    "        \n",
    "        # Deal with impact first (after TC, only)\n",
    "        idx=np.logical_and(yr==locs[ll,1],np.logical_and(grid_times>=tc_jd[ll]+1,grid_times<=tc_jd[ll]+ndays))\n",
    "        n=np.sum(idx)\n",
    "        \n",
    "        if n <ndays: # then we need to \"wrap\" around to the BEGINNING of the NEXT year\n",
    "            idx=np.logical_or(idx,np.logical_and(yr==locs[ll,1]+1,grid_times<=ndays-n))\n",
    "            \n",
    "        # Store impact \n",
    "        impact[v][ll,:]=anoms[v][idx,ll]\n",
    "       \n",
    "        # Now repeat -- but for profile (includes days before -- little more fiddly to set index)\n",
    "        idx_before=np.logical_and(yr==locs[ll,1],np.logical_and(grid_times<=tc_jd[ll],grid_times>=tc_jd[ll]-nbefore))\n",
    "        n=np.sum(idx_before)\n",
    "        if n <(nbefore+1): #then we need to \"wrap\" around to the END of the PREVIOUS year\n",
    "            yr_prev=locs[ll,1]-1\n",
    "            nshort=nbefore-n\n",
    "            if calendar.isleap(np.int(yr_prev)): dmax=366\n",
    "            else: dmax=365  \n",
    "            idx_before=np.logical_or(idx_before,np.logical_and(yr==locs[ll,1]-1,grid_times>=dmax-nshort))\n",
    "        idx=np.logical_or(idx,idx_before)\n",
    "        \n",
    "        # Store profile\n",
    "        profile[v][ll,:]=anoms[v][idx,ll]\n",
    "        nim+=1\n",
    "        \n",
    "    # Truncate as necessary\n",
    "    impact[v]=impact[v][:nim,:]\n",
    "    profile[v]=profile[v][:nim,:]"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## [3] Exploring the results "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 148,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 432x288 with 4 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "# Now we can summarise/explore the profile -- assessing the mean anomaly -- and detecting whether it is significantly \n",
    "# different from zero \n",
    "mu_anom={}\n",
    "x=range(-30,31)\n",
    "fig,ax=plt.subplots(2,2); count=0\n",
    "for v in vs: \n",
    "    a=ax.flat[count]\n",
    "    mu_anom[v]=np.mean(profile[v],axis=0)\n",
    "    sem=np.std(profile[v],axis=0)/np.sqrt(len(profile[v]))\n",
    "    a.plot(x,mu_anom[v])\n",
    "    a.fill_between(x,mu_anom[v]-1.96*sem,mu_anom[v]+1.96*sem,alpha=0.2)\n",
    "    a.axhline(0,color=\"red\"); a.set_title(v); count+=1\n",
    "plt.tight_layout()\n",
    "\n",
    "# Saving the results: (pickling; take anoms and full profile)\n",
    "pickle.dump(anoms,open(\"TC_Met_Anomaly_Time_Series.p\",\"wb\"))\n",
    "pickle.dump(profile,open(\"TC_Met_Anomaly_Profile.p\",\"wb\"))\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# [Extra] Repeat the extraction, but for the hi, only. \n",
    "Here we repeat the same idea as in [2], but we read a different NetCDF4 dataset, which extends up to 2017. And we \n",
    "Also take the raw time series -- not the anomalies. "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 172,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "37837458.02521658\n"
     ]
    }
   ],
   "source": [
    "# Change input data#\n",
    "metfile=\"/media/gytm3/WD12TB/TropicalCyclones/TC-DeadlyHeat/Data/TC_hi_only.nc\"\n",
    "meto=Dataset(metfile,\"r\")\n",
    "\n",
    "# get ntime and nlocs\n",
    "ntime,nlocs=meto.variables[vs[0]].shape\n",
    "time=meto.variables[\"time\"]\n",
    "yr,mon,day,hr,dt=\\\n",
    "GF.conTimes(time_str=time.units,calendar=time.calendar,\\\n",
    "                    times=time[:],safe=False)\n",
    "grid_times=np.floor(tc.dt2decday(dt)) # decimal time (WFDEI)\n",
    "pop=meto.variables[\"population\"][:]\n",
    "\n",
    "# Otherwise, we basically copy/paste/edit the code from [3]...\n",
    "\n",
    "# Make substitutions for clearer code\n",
    "tc_yy=locs[:,1]\n",
    "tc_jd=np.floor(locs[:,2])\n",
    "nim=0\n",
    "hi_out=np.zeros((nlocs,ndays))\n",
    "hi_out_pop=np.zeros((nlocs,3))\n",
    "for ll in range(nlocs):\n",
    "                 \n",
    "        if tc_yy[ll]==2017 and tc_jd[ll]>(365-ndays): continue # because ndays window will trigger an \n",
    "            # out-of-bounds error\n",
    "        \n",
    "        # Deal with impact first (after TC, only)\n",
    "        idx=np.logical_and(yr==locs[ll,1],np.logical_and(grid_times>=tc_jd[ll]+1,grid_times<=tc_jd[ll]+ndays))\n",
    "        n=np.sum(idx)\n",
    "        \n",
    "        if n <ndays: # then we need to \"wrap\" around to the BEGINNING of the NEXT year\n",
    "            idx=np.logical_or(idx,np.logical_and(yr==locs[ll,1]+1,grid_times<=ndays-n))\n",
    "            \n",
    "        # Store raw hi series\n",
    "        hi_out[ll,:]=meto.variables[\"hi\"][idx,ll]\n",
    "        \n",
    "        # Store max hi and population for thos loc\n",
    "        hi_out_pop[ll,0]=ids[ll]; hi_out_pop[ll,1]=np.max(hi_out[ll,:]); hi_out_pop[ll,2]=pop[ll]\n",
    "        nim+=1      \n",
    "              \n",
    "# Truncate as necessary\n",
    "hi_out=hi_out[:nim,:]\n",
    "hi_out_pop=hi_out_pop[:nim,:]\n",
    "\n",
    "# Write the pop file out to a text file (human-readable); the hi_out is post-processed, below.\n",
    "np.savetxt(\"/media/gytm3/WD12TB/TropicalCyclones/TC-DeadlyHeat/Data/HI_pop_out.txt\",hi_out_pop,\n",
    "          fmt=\"%.2f\",header=\"ID\\tHI\\tPop\",comments=\"\")\n",
    "\n",
    "# Now we need to find the maximum across each TC at each time step (using the ID attached to each of the locations)\n",
    "# If we exceed 40.6, we should save the HI series (pickle) using the ID as the filename. Otherwise we produce an \n",
    "# n x 2 array: [id | max hi]\n",
    "ids=meto.variables[\"id\"]\n",
    "uids=np.unique(ids)\n",
    "tc_hi_record=np.zeros((len(uids),2))\n",
    "count=0\n",
    "for ii in uids:\n",
    "    idx=ids==ii\n",
    "    tc_hi_record[count,0]=ii; tc_hi_record[count,1]=np.max(hi_out[idx,:])\n",
    "    if tc_hi_record[count,1]>40.6: \n",
    "        pickle.dump(np.max(hi_out[idx,:],axis=0),open(\"%.0f_max_hi.p\",\"wb\"))\n",
    "    count+=1\n",
    "    \n",
    "# We also need to think about the population at this stage -- should be in the netCDF file...\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 162,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 2",
   "language": "python",
   "name": "python2"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 2
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython2",
   "version": "2.7.15"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
