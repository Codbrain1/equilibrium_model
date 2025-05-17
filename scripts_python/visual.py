from matplotlib.figure import figaspect
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
import numpy as np
from decimal import Decimal, getcontext
getcontext().prec = 30
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
    """Форматирует числа в читаемый вид с фиксированной запятой"""
    if abs(value) < 1e-4 or abs(value) >= 1e6:
        # Для очень малых/больших чисел используем научную нотацию без "e"
        exponent = int(np.floor(np.log10(abs(value))))
        coeff = value / 10**exponent
        return fr"${coeff:.2f} \times 10^{{{exponent}}}$"
    else:
        # Для средних значений используем фиксированную запятую
        return f"{value:.6f}"  # 6 знаков после запятой
    
def format_high_precision(value):
    """Форматирует числа с высокой точностью (9 знаков после запятой)"""
    if value == 0:
        return "0"
    
    # Для значений в диапазоне 1e-9 до 1e9 используем фиксированную запятую
    if 1e-9 <= abs(value) < 1e9:
        return f"{value:.9f}"
    
    # Для очень больших/малых значений используем научную нотацию
    exponent = int(np.floor(np.log10(abs(value))))
    coeff = value / 10**exponent
    return fr"${coeff:.9f} \times 10^{{{exponent}}}$"

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
                t.append(Decimal(parts[0]))
                E.append(Decimal(parts[1]))
                P.append(Decimal(parts[2]))
                L.append(Decimal(parts[3]))
                E_k.append(Decimal(parts[4]))
                E_p.append(Decimal(parts[5]))

    plot_params = {'figsize': (14, 8), 'dpi': 300}
    
    # 1. Полная энергия
    fig, ax = plt.subplots(**plot_params)
    ax.plot(t, E, color='tab:blue', linewidth=2)
    ax.set_xlabel("Время, t")
    ax.set_ylabel("Энергия, E")
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, _: format_scientific_notation(x)))
    ax.grid(True)
    save_plot(fig, "energy.jpg")

    # 2. График относительного изменения энергии
    dEE_0 = [abs(Ei-E[0]+Decimal(1e-31))/abs(E[0]+Decimal(1e-31)) for Ei in E]
    fig, ax = plt.subplots(**plot_params)
    ax.semilogy(t, dEE_0, color='tab:red', linewidth=2)
    ax.set_xlabel("Время, t")
    ax.set_ylabel(r"$\dfrac{\Delta E}{E_0}$")
    ax.grid(True)
    save_plot(fig, "energy_dE_E_0.jpg")

    # 3. График импульса
    fig, ax = plt.subplots(**plot_params)
    ax.plot(t, P, label="Импульс", color='tab:green', linewidth=2)
    ax.set_xlabel("Время")
    ax.set_ylabel("Импульс, P")
    ax.grid(True)
    save_plot(fig, "impulse.jpg")

     
    # 4. График момента импульса
    fig, ax = plt.subplots(**plot_params)
    ax.plot(t, L, label="Момент импульса", color='tab:purple', linewidth=2)
    ax.set_xlabel("Время")
    ax.set_ylabel("Момент, L")
    ax.grid(True)
    save_plot(fig, "moment_impulse.jpg")
    
    # 2. График относительного изменения момента импульса
    dLL_0 = [abs(Li-L[0]+Decimal(1e-31))/(abs(L[0])+Decimal(1e-31)) for Li in L] 
    fig2, ax2 = plt.subplots(**plot_params)
    ax2.semilogy(t, dLL_0, color='tab:red', linewidth=2)
    ax2.set_xlabel("Время")
    ax2.set_ylabel(r"$\dfrac{\Delta L}{L_0}$")
    ax2.grid(True)
    save_plot(fig2, "Momen_impulse_dL_L_0.jpg")
    
    # 5. График кинетической энергии
    fig, ax = plt.subplots(**plot_params)
    ax.plot(t, E_k, label="Кинетическая энергия", color='tab:orange', linewidth=2)
    ax.set_xlabel("Время")
    ax.set_ylabel(r"$E_k$")
    ax.grid(True)
    save_plot(fig, "energy_k.jpg")

    # 6. График потенциальной энергии
    fig, ax = plt.subplots(**plot_params)
    ax.plot(t, E_p, label="Потенциальная энергия", color='tab:brown', linewidth=2)
    ax.set_xlabel("Время")
    ax.set_ylabel(r"$E_p$")
    ax.grid(True)
    save_plot(fig, "energy_p.jpg")
    
def visual_dependens_dt_time():
    configure_plot_settings()
    
    t, dt = [], []
    with open(r"../результат_моделирования/time.txt", 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                t.append(Decimal(parts[0]))
                dt.append(Decimal(parts[1]))
         
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.plot(t, dt, color='tab:blue', linewidth=2)
    ax.set_xlabel("Время работы программы (сек)")
    ax.set_ylabel(r"$\Delta t$")
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, _: format_scientific_notation(x)))
    ax.grid(True)
    save_plot(fig, "time_dependence.jpg")
    
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
            t = Decimal(line.strip())
            for i in range(N):
                parts = f.readline().strip().split()
                if len(parts) >= 2:
                    X[i].append(Decimal(parts[0]))
                    Y[i].append(Decimal(parts[1]))

            fig, ax = plt.subplots(figsize=(12, 10))
            for i in range(N):
                ax.plot(X[i], Y[i], linewidth=1)
                if X[i] and Y[i]:
                    ax.scatter([X[i][-1]], [Y[i][-1]], s=100, c='blue')
        
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(r"$y$")
            ax.grid(True)
            name= "trajectories"+str(ind)+".jpg"
            save_plot(fig, name,consshow=False)
def print_to_cadrs():
    configure_plot_settings()
    count=0
    with open(r"../результат_моделирования/positions.txt", 'r') as f:
        N = int(f.readline())
        for line in f:
            X = [[] for _ in range(N)]
            Y = [[] for _ in range(N)]
            t = Decimal(line.strip())
            for i in range(N):
                parts = f.readline().strip().split()
                if len(parts) >= 2:
                    X[i].append(Decimal(parts[0]))
                    Y[i].append(Decimal(parts[1]))
            fig, ax = plt.subplots(figsize=(12, 10))
            # Рисуем траектории
            for i in range(N):
                ax.plot(X[i], Y[i], linewidth=1)
                if X[i] and Y[i]:
                    ax.scatter([X[i][-1]], [Y[i][-1]], s=2, c='blue')
            
            # Устанавливаем пределы осей (правильный синтаксис)
            ax.set_xlim(-2, 2)  # left=-15, right=15
            ax.set_ylim(-2, 2)  # bottom=-15, top=15
            
            # Исправленная подпись оси Y (была опечатка "sy")
            ax.set_xlabel(r"$x$")
            ax.set_ylabel(r"$y$")
            
            ax.grid(True)
            save_plot(fig, "trajectories"+str(count)+".jpg",consshow=False)
            count+=1
            X.clear()
            Y.clear()
             
def visual_traectories():
    """Визуализация траекторий"""
    configure_plot_settings()
    
    with open(r"../результат_моделирования/positions.txt", 'r') as f:
        N = int(f.readline())
        X = [[] for _ in range(N)]
        Y = [[] for _ in range(N)]
        
        for line in f:
            t = Decimal(line.strip())
            for i in range(N):
                parts = f.readline().strip().split()
                if len(parts) >= 2:
                    X[i].append(Decimal(parts[0]))
                    Y[i].append(Decimal(parts[1]))

    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Рисуем траектории
    for i in range(N):
        ax.plot(X[i], Y[i], linewidth=1)
        if X[i] and Y[i]:
            ax.scatter([X[i][-1]], [Y[i][-1]], s=200, c='blue')
    
    # Устанавливаем пределы осей (правильный синтаксис)
    # ax.set_xlim(-1, 1)  # left=-15, right=15
    # ax.set_ylim(-1, 1)  # bottom=-15, top=15
    # Исправленная подпись оси Y (была опечатка "sy")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.axis("equal")
    ax.grid(True)
    save_plot(fig, "trajectories.jpg")
def visual_centr_mass():
    """Визуализация центра масс"""
    configure_plot_settings()
    
    t, X, Y = [], [], []
    with open(r"../результат_моделирования/centr_mass.txt", 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                t.append(Decimal(parts[0]))
                X.append(Decimal(parts[1]))
                Y.append(Decimal(parts[2]))

    fig, ax = plt.subplots(figsize=(12, 10))
    ax.plot(X, Y, linewidth=1)
    ax.scatter(X[0], Y[0], s=20, c='red')
    ax.scatter(X[-1], Y[-1], s=20, c='green')
    
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.grid(True)
    save_plot(fig, "center_of_mass.jpg")
                
if __name__ == "__main__":
    visual_conversation_laws()
    #visual_traectories()
    visual_dependens_dt_time()
    #visual_centr_mass()
    print_to_cadrs()