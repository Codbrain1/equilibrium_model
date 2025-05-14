from matplotlib.figure import figaspect
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
import numpy as np

def set_scientific_fontsize(ax, fontsize=30):
    """Устанавливает размер шрифта для научной нотации на осях"""
    ax.yaxis.get_offset_text().set_fontsize(fontsize)
    ax.xaxis.get_offset_text().set_fontsize(fontsize)
    
def save_plot(fig, filename, show=False,consshow=True):
    """Унифицированная функция сохранения графиков"""
    path = rf"../результат_моделирования/{filename}"
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    if(consshow==True):
        print(f"График сохранен: {path}")
    if show:
        plt.show()
    plt.close(fig)

def format_scientific_notation(value):
    """Форматирует числа в читаемый вид без 'e' нотации"""
    if value == 0:
        return "0"
    exponent = int(np.floor(np.log10(abs(value))))
    coeff = value / 10**exponent
    if abs(exponent) > 2:
        return fr"${coeff:.2f} \times 10^{{{exponent}}}$"
    return f"{value:.4f}"

def configure_plot_settings():
    """Настройки графиков для всех функций"""
    plt.style.use('default')
    mpl.rcParams.update({
        'axes.titlesize': 32,
        'axes.labelsize': 30,
        'xtick.labelsize': 30,
        'ytick.labelsize': 30,
        'legend.fontsize': 30,
        'figure.titlesize': 30,
        'font.family': 'serif',
        'mathtext.default': 'regular',
        'grid.alpha': 0.2,
        'grid.linestyle': '--',
        'grid.color': 'gray',
        'figure.facecolor': 'white',
        'axes.facecolor': 'white',
        'savefig.facecolor': 'white'
    })

def visual_conversation_laws():
    configure_plot_settings()
    E = []
    P = []
    L = []
    t = []
    E_k = []
    E_p = [] 
    
    #with open(r"D:\научка\результат_моделирования\measurements.txt", 'r') as stream:
    with open(r"../результат_моделирования/measurements.txt", 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 6:
                t.append(float(parts[0]))
                E.append(float(parts[1]))
                P.append(float(parts[2]))
                L.append(float(parts[3]))
                E_k.append(float(parts[4]))
                E_p.append(float(parts[5]))

    plot_params = {'figsize': (14, 8), 'dpi': 300}
    
    # 1. Полная энергия
    fig, ax = plt.subplots(**plot_params)
    ax.plot(t, E, color='tab:blue', linewidth=2)
    ax.set_xlabel("Время, t")
    ax.set_ylabel("Энергия, E")
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, _: format_scientific_notation(x)))
    ax.grid(True)
    save_plot(fig, "energy.jpeg")

    # 2. График относительного изменения энергии
    dEE_0 = [abs(Ei-E[0]+1e-30)/abs(E[0]+1e-30) for Ei in E]
    fig, ax = plt.subplots(**plot_params)
    ax.semilogy(t, dEE_0, color='tab:red', linewidth=2)
    ax.set_xlabel("Время, t")
    ax.set_ylabel(r"$\frac{\Delta E}{E_0}$")
    ax.grid(True)
    save_plot(fig, "energy_dE_E_0.jpeg")

    # 3. График импульса
    fig, ax = plt.subplots(**plot_params)
    ax.plot(t, P, label="Импульс", color='tab:green', linewidth=2)
    ax.set_xlabel("Время")
    ax.set_ylabel("Импульс, P")
    ax.grid(True)
    save_plot(fig, "impulse.jpeg")

     
    # 4. График момента импульса
    fig, ax = plt.subplots(**plot_params)
    ax.plot(t, L, label="Момент импульса", color='tab:purple', linewidth=2)
    ax.set_xlabel("Время")
    ax.set_ylabel("Момент, L")
    ax.grid(True)
    save_plot(fig, "moment_impulse.jpeg")
    
    # 2. График относительного изменения момента импульса
    dLL_0 = [abs(Li-L[0]+1e-30)/(abs(L[0]) + 1e-30) for Li in L] 
    fig2, ax2 = plt.subplots(**plot_params)
    ax2.semilogy(t, dLL_0, color='tab:red', linewidth=2)
    ax2.set_xlabel("Время")
    ax2.set_ylabel(r"$\frac{\Delta L}{L_0}$")
    ax2.grid(True)
    save_plot(fig2, "Momen_impulse_dL_L_0.jpeg")
    
    # 5. График кинетической энергии
    fig, ax = plt.subplots(**plot_params)
    ax.plot(t, E_k, label="Кинетическая энергия", color='tab:orange', linewidth=2)
    ax.set_xlabel("Время")
    ax.set_ylabel(r"$E_k$")
    ax.grid(True)
    save_plot(fig, "energy_k.jpeg")

    # 6. График потенциальной энергии
    fig, ax = plt.subplots(**plot_params)
    ax.plot(t, E_p, label="Потенциальная энергия", color='tab:brown', linewidth=2)
    ax.set_xlabel("Время")
    ax.set_ylabel(r"$E_p$")
    ax.grid(True)
    save_plot(fig, "energy_p.jpeg")
    
def visual_dependens_dt_time():
    configure_plot_settings()
    
    t, dt = [], []
    with open(r"../результат_моделирования/time.txt", 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                t.append(float(parts[0]))
                dt.append(float(parts[1]))
         
        # Настройка глобального стиля

    fig, ax = plt.subplots(figsize=(14, 8))
    ax.plot(t, dt, color='tab:blue', linewidth=2)
    ax.set_xlabel("Время работы программы (сек)")
    ax.set_ylabel(r"$\Delta t$")
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, _: format_scientific_notation(x)))
    ax.grid(True)
    save_plot(fig, "time_dependence.jpeg")
    
def print_to_traectories_cadrs():
    """Визуализация траекторий по кадрам"""
    configure_plot_settings()
    
    with open(r"../результат_моделирования/positions.txt", 'r') as f:
        N = int(f.readline())
        X = [[] for _ in range(N)]
        Y = [[] for _ in range(N)]
        ind=0
        for line in f:
            ind=ind+1
            t = float(line.strip())
            for i in range(N):
                parts = f.readline().strip().split()
                if len(parts) >= 2:
                    X[i].append(float(parts[0]))
                    Y[i].append(float(parts[1]))

            fig, ax = plt.subplots(figsize=(12, 10))
            for i in range(N):
                ax.plot(X[i], Y[i], linewidth=1)
                if X[i] and Y[i]:
                    ax.scatter([X[i][-1]], [Y[i][-1]], s=100, c='blue')
        
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(r"$y$")
            ax.grid(True)
            name= "trajectories"+str(ind)+".jpeg"
            save_plot(fig, name,consshow=False)
def print_to_cadrs():
    count=0 
    plt.style.use('seaborn-v0_8')
    mpl.rcParams.update({
        'axes.titlesize': 32,
        'axes.labelsize': 30,
        'xtick.labelsize': 30,
        'ytick.labelsize': 30,
        'legend.fontsize': 30,
        'figure.titlesize': 30,
        'font.family': 'serif',
        'grid.alpha': 0.2,  # Более светлая сетка
        'grid.linestyle': '--',
        'grid.color': 'gray'  # Серый цвет для сетки
    })
    # Чтение данных из файла
    #with open(r"D:\научка\результат_моделирования\positions.txt", 'r') as stream:
    with open(r"../результат_моделирования/positions.txt", 'r') as stream: 
        N = int(stream.readline())
        while True:
            time = 0
            X = []
            Y = []
            line = stream.readline()
            if not line:
                break
            t = float(line.strip())
            time = t
            
            for i in range(N):
                parts = stream.readline().strip().split()
                if len(parts) != 6:
                    raise ValueError(f"Ожидалось 6 значений для объекта {i}, получено {len(parts)}")
                
                x = float(parts[0])
                y = float(parts[1])
                X.append(x)
                Y.append(y)
                
            div=80
            
            if count%div==0:
                # Построение графика
                fig, ax = plt.subplots(figsize=(12, 10))
                plt.axis('equal')
                #ax.set_aspect('equal', adjustable='box')
                
                # Начальная точка (золотая)
                #ax.scatter(X1, Y1, s=100, c='gold', linewidths=5)
    
                # Конечные точки (синие)
                for i in range(0, N):
                    if X[i] and Y[i]:
                        ax.scatter(X[i], Y[i], s=2, c='blue', linewidths=1.5)
                # ax.set_xlim(-5, 5)  
                # ax.set_ylim(-5, 5)
                fig.suptitle(str(time))        
                # Настройка научной нотации
                ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
                set_scientific_fontsize(ax, 30)
                
                ax.set_xlabel("x", fontsize=30)
                ax.set_ylabel("y", fontsize=30)
                ax.grid(True)
                plt.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.9)
                #output_path = r"D:\научка\результат_моделирования\traectories"+str(count)+".png"
                output_path = f"../../результат_моделирования/traectories{str(count)}.png"
                plt.savefig(output_path)
                plt.close(fig)
                #print(f"График сохранен: {output_path}")
            count+=1
             
def visual_traectories():
    """Визуализация траекторий"""
    configure_plot_settings()
    
    with open(r"../результат_моделирования/positions.txt", 'r') as f:
        N = int(f.readline())
        X = [[] for _ in range(N)]
        Y = [[] for _ in range(N)]
        
        for line in f:
            t = float(line.strip())
            for i in range(N):
                parts = f.readline().strip().split()
                if len(parts) >= 2:
                    X[i].append(float(parts[0]))
                    Y[i].append(float(parts[1]))

    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Рисуем траектории
    for i in range(N):
        ax.plot(X[i], Y[i], linewidth=1)
        if X[i] and Y[i]:
            ax.scatter([X[i][-1]], [Y[i][-1]], s=100, c='blue')
    
    # Устанавливаем пределы осей (правильный синтаксис)
    ax.set_xlim(-15, 15)  # left=-15, right=15
    ax.set_ylim(-15, 15)  # bottom=-15, top=15
    
    # Исправленная подпись оси Y (была опечатка "sy")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    
    ax.grid(True)
    save_plot(fig, "trajectories.jpeg")
def visual_centr_mass():
    """Визуализация центра масс"""
    configure_plot_settings()
    
    t, X, Y = [], [], []
    with open(r"../результат_моделирования/centr_mass.txt", 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                t.append(float(parts[0]))
                X.append(float(parts[1]))
                Y.append(float(parts[2]))

    fig, ax = plt.subplots(figsize=(12, 10))
    ax.plot(X, Y, linewidth=1)
    ax.scatter(X[0], Y[0], s=20, c='red')
    ax.scatter(X[-1], Y[-1], s=20, c='green')
    
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.grid(True)
    save_plot(fig, "center_of_mass.jpeg")
                
if __name__ == "__main__":
    visual_conversation_laws()
    visual_traectories()
    visual_dependens_dt_time()
    visual_centr_mass()
    #print_to_traectories_cadrs()