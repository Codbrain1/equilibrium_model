import matplotlib.pyplot as plt
import os

def set_scientific_fontsize(ax, fontsize=30):
    """Устанавливает размер шрифта для научной нотации на осях"""
    ax.yaxis.get_offset_text().set_fontsize(fontsize)
    ax.xaxis.get_offset_text().set_fontsize(fontsize)

def visual_conversation_laws():
    E = []
    P = []
    M = []
    t = []
    E_k = []
    E_p = []
    #R = []
    
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
            #R.append(_R)
    
    # Настройка стиля для всех графиков
    plt.rcParams.update({
        'axes.titlesize': 30,
        'axes.labelsize': 30,
        'xtick.labelsize': 15,
        'ytick.labelsize': 15,
        'legend.fontsize': 25
    })
    
    # 1. График полной энергии
    fig1, ax1 = plt.subplots(figsize=(12, 8))
    ax1.plot(t, E, label="полная энергия")
    ax1.set_xlabel("t", fontsize=30)
    ax1.set_ylabel("E", fontsize=30)
    ax1.ticklabel_format(axis='y', style='plain', scilimits=(0,0))
    ax1.grid(True)
    ax1.legend(loc='lower right', bbox_to_anchor=(1, 0), ncol=1, framealpha=1)
    plt.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.9)   
    output_path = r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\energy.png"
    plt.savefig(output_path)
    plt.show()
    plt.close(fig1)
    print(f"График сохранен: {output_path}")

    # 2. График импульса
    fig2, ax2 = plt.subplots(figsize=(12, 8))
    ax2.plot(t, P, label="импульс")
    ax2.set_xlabel("t", fontsize=30)
    ax2.set_ylabel("P", fontsize=30)
    ax2.ticklabel_format(axis='y', style='plain', scilimits=(0,0))
    ax2.grid(True)
    ax2.legend(loc='lower right', bbox_to_anchor=(1, 0), ncol=1, framealpha=1)
    plt.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.9)
    output_path = r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\impulse.png"
    plt.savefig(output_path)
    plt.close(fig2)
    print(f"График сохранен: {output_path}")
    
    # 3. График момента импульса
    fig3, ax3 = plt.subplots(figsize=(12, 8))
    ax3.plot(t, M, label="момент импульса")
    ax3.set_xlabel("t", fontsize=30)
    ax3.set_ylabel("M", fontsize=30)
    ax3.ticklabel_format(axis='y', style='plain', scilimits=(0,0))
    ax3.grid(True)
    ax3.legend(loc='upper right', bbox_to_anchor=(1, 1), ncol=1, framealpha=1)
    plt.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.9)
    output_path = r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\moment_impulse.png"
    plt.savefig(output_path)
    plt.show()
    plt.close(fig3)
    print(f"График сохранен: {output_path}")

    # 4. График кинетической энергии
    fig4, ax4 = plt.subplots(figsize=(12, 8))
    ax4.plot(t, E_k, label="кинетическая энергия")
    ax4.set_xlabel("t", fontsize=30)
    ax4.set_ylabel("$E_k$", fontsize=30)
    ax4.ticklabel_format(axis='y', style='plain', scilimits=(0,0))
    ax4.grid(True)
    ax4.legend(loc='lower right', bbox_to_anchor=(1, 0), ncol=1, framealpha=1)
    plt.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.9)
    output_path = r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\energy_k.png"
    plt.savefig(output_path)
    plt.close(fig4)
    print(f"График сохранен: {output_path}")

    # 5. График потенциальной энергии
    fig5, ax5 = plt.subplots(figsize=(12, 8))
    ax5.plot(t, E_p, label="потенциальная энергия")
    ax5.set_xlabel("t", fontsize=30)
    ax5.set_ylabel("$E_p$", fontsize=30)
    ax5.ticklabel_format(axis='y', style='plain', scilimits=(0,0))
  
    ax5.grid(True)
    ax5.legend(loc='lower right', bbox_to_anchor=(1, 0), ncol=1, framealpha=1)
    plt.subplots_adjust(left=0.25, right=0.9, bottom=0.15, top=0.9)
    output_path = r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\energy_p.png"
    plt.savefig(output_path)
    plt.close(fig5)
    print(f"График сохранен: {output_path}")

    # 6. График расстояния
    # fig6, ax6 = plt.subplots(figsize=(12, 8))
    # ax6.plot(t, R, label="расстояние от Солнца до Земли")
    # ax6.set_xlabel("t", fontsize=30)
    # ax6.set_ylabel("r", fontsize=30)
    # ax6.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
    # set_scientific_fontsize(ax6, 30)
    # ax6.grid(True)
    # ax6.legend(loc='lower right', bbox_to_anchor=(1, 0), ncol=1, framealpha=1)
    # plt.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.9)
    # output_path = r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\radius_vector.png"
    # plt.savefig(output_path)
    # plt.close(fig6)
    # print(f"График сохранен: {output_path}")
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
                        ax.scatter([X[i][-1]], [Y[i][-1]], s=100, c='green', linewidths=5)
                
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
                ax.scatter([X[i][-1]], [Y[i][-1]], s=100, c='green', linewidths=5)
                
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
    #visual_traectories()
    #print_to_cadrs()