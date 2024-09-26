import matplotlib.pyplot as plt
import matplotlib.animation as animation
import cartopy as cy
import numpy as np
from cmcrameri import cm
import xarray as xr


def findindex(alat,alon,plat,plon):
    #finding identical location of pos plat, plon in array alat[],alon[]
    abslat = np.abs(alat-plat)
    abslon = np.abs(alon-plon)
    c = np.maximum(abslon,abslat)
    y, x = np.where(c == np.min(c))
    #print(alats[x,y],alon[x,y])
    y=int(y[0])
    x=int(x[0])
    return (y,x)


day='26'; month='10'; year='2014'; HH='12'

url='https://thredds.met.no/thredds/dodsC/aromemetcoopstarc/'+str(year)+'/'+str(month)+'/'+str(day)+'/AROME_MetCoOp_'+str(HH)+'_fp.nc_'+str(year)+str(month)+str(day)
arome = xr.open_dataset(url)

url='https://thredds.met.no/thredds/dodsC/metusers/maltem/GEO4902_2020/Arctic.ECMWF_extracted_'+str(year)+str(month)+str(day)+'T'+str(HH)+'Z.nc'
ecifs =  xr.open_dataset(url)


lat0=60.8608; lon0=7.1118 # Flåm in Western Norway

[yloc,xloc] = findindex(arome.latitude,arome.longitude,lat0,lon0) # from regional model AROME MetCoOp



def create_animation():
    levelsPP = range(0, 140, 10)
    projection = cy.crs.LambertConformal(
        central_longitude=arome.projection_lambert.longitude_of_central_meridian,
        central_latitude=arome.projection_lambert.latitude_of_projection_origin,
        standard_parallels=arome.projection_lambert.standard_parallel
    )

    # Adjust figure size to reduce vertical space
    fig, axsm = plt.subplots(1, 2, subplot_kw={'projection': projection}, figsize=(15, 8), sharex=True, sharey=True)
    
    # Adjust subplot layout to reduce vertical whitespace
    plt.subplots_adjust(top=0.9, bottom=0.15)

    for ax in axsm.flatten():
        ax.set_extent([lon0-2, lon0+1.5, lat0-1.1, lat0+1.1], cy.crs.PlateCarree())
        ax.add_feature(cy.feature.COASTLINE, alpha=0.5)
        ax.add_feature(cy.feature.BORDERS, alpha=0.5)
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
        ax.plot(lon0, lat0, 'ro', markersize=5, transform=cy.crs.PlateCarree())
        ax.text(lon0 + 0.1, lat0, 'Flåm', fontsize=14, color='red', transform=cy.crs.PlateCarree())

    # Adjust colorbar position
    cbaxes = fig.add_axes([0.1, 0.08, 0.8, 0.03])
    
    def animate(timestep):
        axsm[0].clear()
        axsm[1].clear()
        
        for ax in axsm:
            ax.set_extent([lon0-2, lon0+1.5, lat0-1.1, lat0+1.1], cy.crs.PlateCarree())
            ax.add_feature(cy.feature.COASTLINE, alpha=0.5)
            ax.add_feature(cy.feature.BORDERS, alpha=0.5)
            gl = ax.gridlines(draw_labels=True)
            gl.top_labels = False
            gl.right_labels = False
            ax.plot(lon0, lat0, 'ro', markersize=5, transform=cy.crs.PlateCarree())
            ax.text(lon0 + 0.1, lat0, 'Flåm', fontsize=14, color='red', transform=cy.crs.PlateCarree())
        
        arome_plot = arome.precipitation_amount_acc.isel(time=timestep, height0=0).plot.contourf(
            ax=axsm[0], transform=projection, cmap=cm.davos_r, extend='max', levels=levelsPP, add_colorbar=False
        )
        axsm[0].set_title(f'AROME MetCoOp {np.datetime_as_string(arome.time[timestep], unit="h")}')

        ecifs_plot = (ecifs.TP.isel(time=timestep) * 1000).plot.contourf(
            ax=axsm[1], transform=cy.crs.PlateCarree(), cmap=cm.davos_r, extend='max', levels=levelsPP, add_colorbar=False
        )
        axsm[1].set_title(f'ECMWF IFS {np.datetime_as_string(arome.time[timestep], unit="h")}')

        if timestep == 0:
            cbar = plt.colorbar(ecifs_plot, cax=cbaxes, orientation='horizontal', extend='max')
            cbar.set_label('Accumulated precipitation [mm]')
        
        return axsm

    anim = animation.FuncAnimation(fig, animate, frames=49, interval=200, blit=False)
    
    # Save the animation
    anim.save('precipitation_animation.gif', writer='pillow', fps=5)
    
    plt.close(fig)

# Call the function to create and save the animation
# create_animation()



def create_temperature_wind_animation():
    Trange = range(-8, 15, 1)
    projection = cy.crs.LambertConformal(
        central_longitude=arome.projection_lambert.longitude_of_central_meridian,
        central_latitude=arome.projection_lambert.latitude_of_projection_origin,
        standard_parallels=arome.projection_lambert.standard_parallel
    )

    # Adjust figure size to reduce vertical space
    fig, axsm = plt.subplots(1, 2, subplot_kw={'projection': projection}, figsize=(15, 8), sharex=True, sharey=True)
    
    # Adjust subplot layout to reduce vertical whitespace
    plt.subplots_adjust(top=0.9, bottom=0.15)

    for ax in axsm.flatten():
        ax.set_extent([lon0-2, lon0+1.5, lat0-1.1, lat0+1.1], cy.crs.PlateCarree())
        ax.add_feature(cy.feature.COASTLINE, alpha=0.5)
        ax.add_feature(cy.feature.BORDERS, alpha=0.5)
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
        ax.plot(lon0, lat0, 'go', markersize=5, transform=cy.crs.PlateCarree())
        ax.text(lon0 + 0.1, lat0, 'Flåm', fontsize=14, color='tab:green', transform=cy.crs.PlateCarree())

    # Adjust colorbar position
    cbaxes = fig.add_axes([0.1, 0.08, 0.8, 0.03])

    def animate(timestep):
        axsm[0].clear()
        axsm[1].clear()
        
        for ax in axsm:
            ax.set_extent([lon0-2, lon0+1.5, lat0-1.1, lat0+1.1], cy.crs.PlateCarree())
            ax.add_feature(cy.feature.COASTLINE, alpha=0.5)
            ax.add_feature(cy.feature.BORDERS, alpha=0.5)
            gl = ax.gridlines(draw_labels=True)
            gl.top_labels = False
            gl.right_labels = False
            ax.plot(lon0, lat0, 'go', markersize=5, transform=cy.crs.PlateCarree())
            ax.text(lon0 + 0.1, lat0, 'Flåm', fontsize=14, color='tab:green', transform=cy.crs.PlateCarree())

        # AROME plot
        arome_temp_C = arome.air_temperature_2m.isel(time=timestep, height1=0) - 273.15
        arome_plot = arome_temp_C.plot.contourf(ax=axsm[0], transform=projection, cmap='inferno', extend='max', levels=Trange, add_colorbar=False)
        axsm[0].set_title(f'AROME MetCoOp {np.datetime_as_string(arome.time[timestep], unit="h")}')

        # ECMWF IFS plot
        ecifs_plot = (ecifs.T2M.isel(time=timestep) - 273.15).plot.contourf(ax=axsm[1], transform=cy.crs.PlateCarree(), cmap='inferno', extend='max', levels=Trange, add_colorbar=False)
        axsm[1].set_title(f'ECMWF IFS {np.datetime_as_string(arome.time[timestep], unit="h")}')

        # Add wind for AROME
        u_arome = arome.x_wind_10m.isel(time=timestep, height2=0).values
        v_arome = arome.y_wind_10m.isel(time=timestep, height2=0).values
        lat_arome = arome.latitude.values
        lon_arome = arome.longitude.values
        step = 4
        axsm[0].quiver(lon_arome[::step, ::step], lat_arome[::step, ::step], u_arome[::step, ::step], v_arome[::step, ::step],
                       transform=cy.crs.PlateCarree(), scale=400, zorder=3, color='white')

        # Add wind for ECMWF IFS
        u_ecifs = ecifs.U10M.isel(time=timestep).values
        v_ecifs = ecifs.V10M.isel(time=timestep).values
        lat_ecifs_1d = ecifs.lat.values
        lon_ecifs_1d = ecifs.lon.values
        lat_ecifs, lon_ecifs = np.meshgrid(lat_ecifs_1d, lon_ecifs_1d, indexing='ij')
        step = 1
        axsm[1].quiver(lon_ecifs[::step, ::step], lat_ecifs[::step, ::step], u_ecifs[::step, ::step], v_ecifs[::step, ::step],
                       transform=cy.crs.PlateCarree(), scale=400, zorder=3, color='white')

        if timestep == 0:
            cbar = plt.colorbar(ecifs_plot, cax=cbaxes, orientation='horizontal', extend='max', label='2m Temperature [°C]')

        return axsm

    anim = animation.FuncAnimation(fig, animate, frames=len(arome.time), interval=200, blit=False)
    
    # Save the animation
    anim.save('temperature_wind_animation8.gif', writer='pillow', fps=5)
    
    plt.close(fig)

# Call the function to create and save the animation
# create_temperature_wind_animation()




def polar_low20200204_animation():
    Prange = range(98000, 100000, 50)
    projection = cy.crs.LambertConformal(
        central_longitude=arome.projection_lambert.longitude_of_central_meridian,
        central_latitude=arome.projection_lambert.latitude_of_projection_origin,
        standard_parallels=arome.projection_lambert.standard_parallel
    )

    # Adjust figure size to reduce vertical space
    fig, axsm = plt.subplots(1, 2, subplot_kw={'projection': projection}, figsize=(15, 8), sharex=True, sharey=True)
    
    # Adjust subplot layout to reduce vertical whitespace
    plt.subplots_adjust(top=0.9, bottom=0.15)

    for ax in axsm.flatten():
        ax.set_extent([-1, 20, 57, 73], cy.crs.PlateCarree())
        ax.add_feature(cy.feature.COASTLINE, alpha=0.5)
        ax.add_feature(cy.feature.BORDERS, alpha=0.5)
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False

    # Adjust colorbar position
    cbaxes = fig.add_axes([0.1, 0.08, 0.8, 0.03])

    def animate(timestep):
        axsm[0].clear()
        axsm[1].clear()
        
        for ax in axsm:
            ax.set_extent([-1, 20, 56, 73], cy.crs.PlateCarree())
            ax.add_feature(cy.feature.COASTLINE, alpha=0.5)
            ax.add_feature(cy.feature.BORDERS, alpha=0.5)
            gl = ax.gridlines(draw_labels=True)
            gl.top_labels = False
            gl.right_labels = False

        # AROME plot
        try:
            arome_temp_C = arome.air_pressure_at_sea_level.isel(time=timestep, ensemble_member=1 ,height_above_msl=0)
        except:
            arome_temp_C = arome.air_pressure_at_sea_level.isel(time=timestep, height_above_msl=0)
        
        arome_plot = arome_temp_C.plot.contourf(ax=axsm[0], transform=projection, cmap=cm.davos, extend='max', levels=Prange, add_colorbar=False)
        axsm[0].set_title(f'AROME MetCoOp {np.datetime_as_string(arome.time[timestep], unit="h")}')

        # ECMWF IFS plot
        ecifs_plot = (ecifs.SP.isel(time=timestep)).plot.contourf(ax=axsm[1], transform=cy.crs.PlateCarree(), cmap=cm.davos, extend='max', levels=Prange, add_colorbar=False)
        axsm[1].set_title(f'ECMWF IFS {np.datetime_as_string(arome.time[timestep], unit="h")}')

        if timestep == 0:
            cbar = plt.colorbar(ecifs_plot, cax=cbaxes, orientation='horizontal', extend='max', label='Surface Pressure [Pa]')

        return axsm

    anim = animation.FuncAnimation(fig, animate, frames=48, interval=200, blit=False)
    
    # Save the animation
    anim.save('polar_low_animation.gif', writer='pillow', fps=5)
    
    plt.close(fig)


day='03'; month='02'; year='2020'; HH='00'

url='https://thredds.met.no/thredds/dodsC/aromemetcoopstarc/'+str(year)+'/'+str(month)+'/'+str(day)+'/AROME_MetCoOp_'+str(HH)+'_fp.nc_'+str(year)+str(month)+str(day)

working_url = 'https://thredds.met.no/thredds/dodsC/mepsoldarchive/2020/02/03/meps_full_2_5km_20200203T00Z.nc'
arome = xr.open_dataset(working_url)

url='https://thredds.met.no/thredds/dodsC/metusers/maltem/GEO4902_2020/Arctic.ECMWF_extracted_'+str(year)+str(month)+str(day)+'T'+str(HH)+'Z.nc'
ecifs =  xr.open_dataset(url)

# Call the function to create and save the animation
polar_low20200204_animation()



def polar_low20200204_wind_animation():
    Prange = range(98000, 100000, 50)
    projection = cy.crs.LambertConformal(
        central_longitude=arome.projection_lambert.longitude_of_central_meridian,
        central_latitude=arome.projection_lambert.latitude_of_projection_origin,
        standard_parallels=arome.projection_lambert.standard_parallel
    )

    # Adjust figure size to reduce vertical space
    fig, axsm = plt.subplots(1, 2, subplot_kw={'projection': projection}, figsize=(15, 8), sharex=True, sharey=True)
    
    # Adjust subplot layout to reduce vertical whitespace
    plt.subplots_adjust(top=0.9, bottom=0.15)

    for ax in axsm.flatten():
        ax.set_extent([-1, 20, 57, 73], cy.crs.PlateCarree())
        ax.add_feature(cy.feature.COASTLINE, alpha=0.5)
        ax.add_feature(cy.feature.BORDERS, alpha=0.5)
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False

    # Adjust colorbar position
    cbaxes = fig.add_axes([0.1, 0.08, 0.8, 0.03])

    def animate(timestep):
        axsm[0].clear()
        axsm[1].clear()
        
        for ax in axsm:
            ax.set_extent([-1, 9, 57.5, 64], cy.crs.PlateCarree())
            ax.add_feature(cy.feature.COASTLINE, alpha=0.5)
            ax.add_feature(cy.feature.BORDERS, alpha=0.5)
            gl = ax.gridlines(draw_labels=True)
            gl.top_labels = False
            gl.right_labels = False

            ax.plot(lon0, lat0, 'ro', markersize=5, transform=cy.crs.PlateCarree())
            ax.text(lon0 + 0.1, lat0, 'Flåm', fontsize=14, color='red', transform=cy.crs.PlateCarree())

        # AROME plot
        try:
            arome_temp_C = arome.air_pressure_at_sea_level.isel(time=timestep, ensemble_member=1 ,height_above_msl=0)
        except:
            arome_temp_C = arome.air_pressure_at_sea_level.isel(time=timestep, height_above_msl=0)
        
        arome_plot = arome_temp_C.plot.contourf(ax=axsm[0], transform=projection, cmap=cm.davos, extend='max', levels=Prange, add_colorbar=False)
        axsm[0].set_title(f'AROME MetCoOp {np.datetime_as_string(arome.time[timestep], unit="h")}')

        # ECMWF IFS plot
        ecifs_plot = (ecifs.SP.isel(time=timestep)).plot.contourf(ax=axsm[1], transform=cy.crs.PlateCarree(), cmap=cm.davos, extend='max', levels=Prange, add_colorbar=False)
        axsm[1].set_title(f'ECMWF IFS {np.datetime_as_string(arome.time[timestep], unit="h")}')

        # Add wind for AROME
        u_arome = arome.x_wind_10m.isel(time=timestep, height2=0).values
        v_arome = arome.y_wind_10m.isel(time=timestep, height2=0).values
        lat_arome = arome.latitude.values
        lon_arome = arome.longitude.values
        step = 4
        axsm[0].quiver(lon_arome[::step, ::step], lat_arome[::step, ::step], u_arome[::step, ::step], v_arome[::step, ::step],
                       transform=cy.crs.PlateCarree(), scale=400, zorder=3, color='white')

        # Add wind for ECMWF IFS
        u_ecifs = ecifs.U10M.isel(time=timestep).values
        v_ecifs = ecifs.V10M.isel(time=timestep).values
        lat_ecifs_1d = ecifs.lat.values
        lon_ecifs_1d = ecifs.lon.values
        lat_ecifs, lon_ecifs = np.meshgrid(lat_ecifs_1d, lon_ecifs_1d, indexing='ij')
        step = 1
        axsm[1].quiver(lon_ecifs[::step, ::step], lat_ecifs[::step, ::step], u_ecifs[::step, ::step], v_ecifs[::step, ::step],
                       transform=cy.crs.PlateCarree(), scale=400, zorder=3, color='white')

        if timestep == 0:
            cbar = plt.colorbar(ecifs_plot, cax=cbaxes, orientation='horizontal', extend='max', label='Surface Pressure [Pa]')

        return axsm

    anim = animation.FuncAnimation(fig, animate, frames=48, interval=200, blit=False)
    
    # Save the animation
    anim.save('polar_low_animation_zoom.gif', writer='pillow', fps=5)
    
    plt.close(fig)


day='03'; month='02'; year='2020'; HH='00'

url='https://thredds.met.no/thredds/dodsC/aromemetcoopstarc/'+str(year)+'/'+str(month)+'/'+str(day)+'/AROME_MetCoOp_'+str(HH)+'_fp.nc_'+str(year)+str(month)+str(day)

working_url = 'https://thredds.met.no/thredds/dodsC/mepsoldarchive/2020/02/03/meps_full_2_5km_20200203T00Z.nc'
arome = xr.open_dataset(working_url)

url='https://thredds.met.no/thredds/dodsC/metusers/maltem/GEO4902_2020/Arctic.ECMWF_extracted_'+str(year)+str(month)+str(day)+'T'+str(HH)+'Z.nc'
ecifs =  xr.open_dataset(url)

# Call the function to create and save the animation
polar_low20200204_wind_animation()