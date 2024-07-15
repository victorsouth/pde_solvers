import pandas as pd
from matplotlib import pyplot as plt
import os 

current_dir = os.path.abspath(os.getcwd())

main_path = os.path.dirname(os.path.dirname(current_dir))

#research_path = '..\\..\\research_out\\research_out\\QuasiStationaryModel\\ShowProfileImpactInQuasiStationaryModel\\'
folders = os.path.join('research_out', 'QSM_models', 'QuasiStationaryModel', 'ShowProfileImpactInQuasiStationaryModel')
case_paths = ['path_full_profile', 'start_end_profile']
file_name = 'output pressure.csv'
coordinates = 'output coordinates.csv'
heights = 'output heights.csv'
pr = 'output profile.csv'
plotsCount = 2




def plot_coordinates_heights(ax):
    profiles = []
    for path in case_paths:
        profile_path = os.path.join(main_path, folders, path, pr)
        prof = pd.read_csv(profile_path, header=None, sep=';').T
        prof = prof[1:-1]
        prof.columns = ['coord', 'heights']
        print(prof)
        profiles.append(prof)
    ax.plot(profiles[0]['coord'], profiles[0]['heights'])
    ax.plot(profiles[1]['coord'], profiles[1]['heights'])
    ax.set_xlabel('Координаты')
    ax.set_ylabel('Высотки')
    ax.grid()
        
def plot_diff_press_at_the_end(ax):
    profiles = []
    
    for path in case_paths:
        research_path = os.path.join(main_path, folders, path, file_name)
        df = pd.read_csv(research_path, header=None, sep=';')
        df[df.columns[0]] = df[df.columns[0]].astype('datetime64[ns]')
        profiles.append(df.iloc[:, :-1])
    
    res = profiles[1].iloc[:,-1].astype('float') - profiles[0].iloc[:,-1].astype('float')
    print(df[df.columns[0]], res)
    
    ax.plot(df[df.columns[0]], res)
    ax.set_xlabel('Время')
    ax.set_ylabel('Разница давления в конце трубопровода')    
    ax.grid()   
    
if __name__ == '__main__':
    axes = [plt.subplot(plotsCount, 1, _ + 1) for _ in range(plotsCount)]
    plot_coordinates_heights(axes[0])   
    plot_diff_press_at_the_end(axes[1])
    plt.show()    