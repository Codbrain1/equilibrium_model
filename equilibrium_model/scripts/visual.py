from matplotlib.figure import figaspect
import matplotlib.pyplot as plt
import os
import matplotlib as mpl

def set_scientific_fontsize(ax, fontsize=30):
    """Устанавливает размер шрифта для научной нотации на осях"""
    ax.yaxis.get_offset_text().set_fontsize(fontsize)
    ax.xaxis.get_offset_text().set_fontsize(fontsize)
    
def save_plot(fig, filename, show=False):
    """Унифицированная функция сохранения графиков"""
    path = rf"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\{filename}"
    fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"График сохранен: {path}")
    if show:
        plt.show()
    plt.close(fig)
def visual_conversation_laws():
    E = []
    P = []
    M = []
    t = []
    E_k = []
    E_p = [] 
    with open(r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\measurements.txt", 'r') as stream:
        while True:
            line = stream.readline()
            if not line:
                break
                
            parts = line.strip().split()
            if len(parts) != 6:
                continue
                
            time = float(parts[0])
            _E = float(parts[1])
            _P = float(parts[2])
            _M = float(parts[3])
            _E_k = float(parts[4])
            _E_p = float(parts[5])
            #_R = float(parts[6])
            
            E.append(_E)
            P.append(_P)
            M.append(_M)
            t.append(time)
            E_k.append(_E_k)
            E_p.append(_E_p)

   # Глобальные настройки стиля
    plt.style.use('seaborn-v0_8')
    mpl.rcParams.update({
        'axes.titlesize': 24,
        'axes.labelsize': 22,
        'xtick.labelsize': 18,
        'ytick.labelsize': 18,
        'legend.fontsize': 20,
        'figure.titlesize': 26,
        'font.family': 'serif',
        'grid.alpha': 0.2,  # Более светлая сетка
        'grid.linestyle': '--',
        'grid.color': 'gray'  # Серый цвет для сетки
    })

    # Общие параметры для всех графиков
    plot_params = {
        'figsize': (14, 8),
        'dpi': 100
    }

    # 1. График полной энергии
    fig1, ax1 = plt.subplots(**plot_params)
    ax1.plot(t, E, label="Полная энергия", color='tab:blue', linewidth=2)
    ax1.set_xlabel("Время, t", fontweight='bold')
    ax1.set_ylabel("Энергия, E", fontweight='bold')
    ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    ax1.grid(True)  # Сетка включается отдельно для каждого axes
    ax1.legend(loc='upper right', framealpha=0.9)
    save_plot(fig1, "energy.png", show=True)

    # 2. График относительного изменения энергии
    dEE_0 = [abs(Ei-E[0])/abs(E[0]) for Ei in E]
    fig2, ax2 = plt.subplots(**plot_params)
    ax2.semilogy(t, dEE_0, label="KDK", color='tab:red', linewidth=2)
    ax2.set_xlabel("Время, t", fontweight='bold')
    ax2.set_ylabel(r"$\frac{\Delta E}{E_0}$", fontweight='bold')
    ax2.grid(True)
    ax2.legend(loc='upper right', framealpha=0.9)
    save_plot(fig2, "energy_dE_E_0.png")

    # 3. График импульса
    fig3, ax3 = plt.subplots(**plot_params)
    ax3.plot(t, P, label="Импульс", color='tab:green', linewidth=2)
    ax3.set_xlabel("Время, t", fontweight='bold')
    ax3.set_ylabel("Импульс, P", fontweight='bold')
    ax3.grid(True)
    ax3.legend(loc='upper right', framealpha=0.9)
    save_plot(fig3, "impulse.png")

    # 4. График момента импульса
    fig4, ax4 = plt.subplots(**plot_params)
    ax4.plot(t, M, label="Момент импульса", color='tab:purple', linewidth=2)
    ax4.set_xlabel("Время, t", fontweight='bold')
    ax4.set_ylabel("Момент, M", fontweight='bold')
    ax4.grid(True)
    ax4.legend(loc='upper right', framealpha=0.9)
    save_plot(fig4, "moment_impulse.png", show=True)

    # 5. График кинетической энергии
    fig5, ax5 = plt.subplots(**plot_params)
    ax5.plot(t, E_k, label="Кинетическая энергия", color='tab:orange', linewidth=2)
    ax5.set_xlabel("Время, t", fontweight='bold')
    ax5.set_ylabel(r"$E_k$", fontweight='bold')
    ax5.grid(True)
    ax5.legend(loc='upper right', framealpha=0.9)
    save_plot(fig5, "energy_k.png")

    # 6. График потенциальной энергии
    fig6, ax6 = plt.subplots(**plot_params)
    ax6.plot(t, E_p, label="Потенциальная энергия", color='tab:brown', linewidth=2)
    ax6.set_xlabel("Время, t", fontweight='bold')
    ax6.set_ylabel(r"$E_p$", fontweight='bold')
    ax6.grid(True)
    ax6.legend(loc='upper right', framealpha=0.9)
    save_plot(fig6, "energy_p.png")
    
def visual_dependens_dt_time():
    time = []
    dt = []
    ind=0
    with open(r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\time.txt",'r') as stream:
     while True:
         line = stream.readline()
         if not line:
             break
         parts = line.strip().split()
         if(len(parts)!=2):
            print("\033[3m\033[33m\033[41m{}\033[0m".format(f"Ошибка неправильная запись в файле time.txt, строка {ind}"))
            continue
         
         _time = float(parts[0])
         _dt = float(parts[1])
         
         time.append(_time)
         dt.append(_dt)
         ind=ind+1
         
        # Настройка глобального стиля
    plt.style.use('seaborn-v0_8')
    mpl.rcParams.update({
        'axes.titlesize': 24,
        'axes.labelsize': 22,
        'xtick.labelsize': 18,
        'ytick.labelsize': 18,
        'legend.fontsize': 20,
        'figure.titlesize': 26,
        'font.family': 'serif',
        'grid.alpha': 0.2,  # Более светлая сетка
        'grid.linestyle': '--',
        'grid.color': 'gray'  # Серый цвет для сетки
    })

    # Создаем график с новыми параметрами
    fig1, ax1 = plt.subplots(figsize=(14, 8))

    # Основной график
    ax1.plot(time, dt, label="Шаг интегрирования", 
             color='tab:blue', linewidth=2.5)

    # Настройки осей
    ax1.set_xlabel("Время, t", fontsize=22, fontweight='bold')
    ax1.set_ylabel("Δt", fontsize=22, fontweight='bold')
    ax1.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

    # Настройка сетки (более светлая)
    ax1.grid(True, linestyle=':', alpha=0.3, color='lightgray')

    # Оптимизированная легенда
    legend = ax1.legend(
        loc='upper right',        # Новое расположение
        frameon=True,
        framealpha=0.9,
        edgecolor='gray',        # Более светлая граница
        facecolor='whitesmoke',  # Светлый фон
        borderpad=1.2,
        shadow=True              # Легкая тень для объема
    )

    # Настройка расположения графика
    plt.tight_layout(pad=2.0)  # Увеличенный отступ

    # Сохранение с улучшенными параметрами
    output_path = r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\time.png"
    plt.savefig(output_path, 
               dpi=350,                  # Улучшенное качество
               bbox_inches='tight', 
               facecolor='white',        # Белый фон при сохранении
               edgecolor='none')

    plt.show()
    plt.close(fig1)
    print(f"График сохранен: {output_path}")
    
def print_to_cadrs():
    count=0
    # Чтение данных из файла
    with open(r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\positions.txt", 'r') as stream:
        N = int(stream.readline())
        time = []
        X = [[] for _ in range(N)]
        Y = [[] for _ in range(N)]
        X1 = [0]
        Y1 = [0]
        
        while True:
            line = stream.readline()
            if not line:
                break
            t = float(line.strip())
            time.append(t)
            
            for i in range(N):
                parts = stream.readline().strip().split()
                if len(parts) != 6:
                    raise ValueError(f"Ожидалось 6 значений для объекта {i}, получено {len(parts)}")
                
                x = float(parts[0])
                y = float(parts[1])
                X[i].append(x)
                Y[i].append(y)
                
            div=100
            
            if count%div==0:
                # Построение графика
                fig, ax = plt.subplots(figsize=(12, 10))
                plt.axis('equal')
                #ax.set_aspect('equal', adjustable='box')
                # Настройки шрифтов
                plt.rcParams.update({
                    'axes.labelsize': 30,
                    'xtick.labelsize': 30,
                    'ytick.labelsize': 30,
                })
                
                
                # Начальная точка (золотая)
                #ax.scatter(X1, Y1, s=100, c='gold', linewidths=5)
    
                # Траектории
                for i in range(N):
                    ax.plot(X[i], Y[i], linewidth=1)
    
                # Конечные точки (зеленые)
                for i in range(0, N):
                    if X[i] and Y[i]:
                        ax.scatter([X[i][-1]], [Y[i][-1]], s=2, c='blue', linewidths=1.5)
                
                # ax.set_xlim(-5, 5)  
                # ax.set_ylim(-5, 5)
                fig.suptitle(str(time[-1]))        
                # Настройка научной нотации
                ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
                set_scientific_fontsize(ax, 30)
                
                ax.set_xlabel("x", fontsize=30)
                ax.set_ylabel("y", fontsize=30)
                ax.grid(True)
                plt.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.9)
                output_path = r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\traectories"+str(count)+".png"
                plt.savefig(output_path)
                plt.close(fig)
                #print(f"График сохранен: {output_path}")
            count+=1
            
def visual_traectories():
    # Чтение данных из файла
    with open(r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\positions.txt", 'r') as stream:
        N = int(stream.readline())
        time = []
        X = [[] for _ in range(N)]
        Y = [[] for _ in range(N)]
        X1 = [0]
        Y1 = [0]
        
        while True:
            line = stream.readline()
            if not line:
                break
            t = float(line.strip())
            time.append(t)
            
            for i in range(N):
                parts = stream.readline().strip().split()
                if len(parts) != 6:
                    raise ValueError(f"Ожидалось 6 значений для объекта {i}, получено {len(parts)}")
                
                x = float(parts[0])
                y = float(parts[1])
                X[i].append(x)
                Y[i].append(y)
                

        # Построение графика
        fig, ax = plt.subplots(figsize=(12, 10))
        plt.axis('equal')
        # Настройки шрифтов
        plt.rcParams.update({
            'axes.labelsize': 30,
            'xtick.labelsize': 30,
            'ytick.labelsize': 30,
        })

            # Начальная точка (золотая)
            #ax.scatter(X1, Y1, s=100, c='gold', linewidths=5)
    
        # Траектории
        for i in range(N):
            ax.plot(X[i], Y[i], linewidth=1)
    
        # Конечные точки (зеленые)
        for i in range(0, N):
            if X[i] and Y[i]:
                ax.scatter([X[i][-1]], [Y[i][-1]], s=1.5, c='blue', linewidths=1.5)
                
        # ax.set_xlim(-5, 5)  
        # ax.set_ylim(-5, 5)
        fig.suptitle(str(time[-1]))        
        # Настройка научной нотации
        ax.ticklabel_format(axis='both', style='sci', scilimits=(0,0))
        set_scientific_fontsize(ax, 30)
                
        ax.set_xlabel("x", fontsize=30)
        ax.set_ylabel("y", fontsize=30)
        ax.grid(True)
        plt.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.9)
        output_path = r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\traectories.png"
        plt.savefig(output_path)
        plt.close(fig)
        print(f"График сохранен: {output_path}")

if __name__ == "__main__":
    visual_conversation_laws()
    visual_traectories()
    # print_to_cadrs()
    visual_dependens_dt_time()