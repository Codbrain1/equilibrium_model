import matplotlib.pyplot as plt
import os
def visual_converation_laws():
    E = []
    P = []
    M = []
    t = []
    E_k = []
    E_p = []
    R = []
    
    with open(r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\measurements.txt", 'r') as stream:
        while True:
            line = stream.readline()
            if not line:
                break
                
            parts = line.strip().split()
            if len(parts) != 7:
                continue
                
            time = float(parts[0])
            _E = float(parts[1])
            _P = float(parts[2])
            _M = float(parts[3])
            _E_k = float(parts[4])
            _E_p = float(parts[5])
            _R = float(parts[6])
            
            E.append(_E)
            P.append(_P)
            M.append(_M)
            t.append(time)
            E_k.append(_E_k)
            E_p.append(_E_p)
            R.append(_R)
    
    plt.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.9, wspace=0.4, hspace=0.4)
    plt.plot(t, E,label="полная энергия")
    plt.xlabel("t", fontsize=20)
    plt.ylabel("E", fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid(True)
    plt.legend(loc='lower right', bbox_to_anchor=(1, 0), 
           ncol=1, fontsize=12, framealpha=1)
    plt.savefig(r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\energy.png")
    
    plt.figure(2)
    plt.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.9, wspace=0.4, hspace=0.4)
    plt.plot(t, P,label="импульс")
    plt.xlabel("t", fontsize=20)
    plt.ylabel("P", fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid(True)
    plt.legend(loc='lower right', bbox_to_anchor=(1, 0), 
           ncol=1, fontsize=12, framealpha=1)
    plt.savefig(r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\impulse.png")
    
    plt.figure(3)
    plt.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.9, wspace=0.4, hspace=0.4)
    plt.plot(t, M,label="момент импульса")
    plt.xlabel("t", fontsize=20)
    plt.ylabel("M", fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid(True)
    plt.legend(loc='lower right', bbox_to_anchor=(1, 0), 
           ncol=1, fontsize=12, framealpha=1)
    plt.savefig(r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\moment_impulse.png")
    
    plt.figure(4)
    plt.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.9, wspace=0.4, hspace=0.4)
    plt.plot(t, E_k,label="кинетическая энергия")
    plt.xlabel("t", fontsize=20)
    plt.ylabel("$E_k$", fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid(True)
    plt.legend(loc='lower right', bbox_to_anchor=(1, 0), 
           ncol=1, fontsize=12, framealpha=1)
    plt.savefig(r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\energy_k.png")
    
    plt.figure(5)
    plt.subplots_adjust(left=0.25, right=0.9, bottom=0.15, top=0.9, wspace=0.4, hspace=0.4)
    plt.plot(t, E_p,label="потенциальная энергия")
    plt.xlabel("t", fontsize=20)
    plt.ylabel("$E_p$", fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid(True)
    plt.legend(loc='lower right', bbox_to_anchor=(1, 0), 
           ncol=1, fontsize=12, framealpha=1)
    plt.savefig(r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\energy_p.png")
    
    plt.figure(6)
    plt.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.9, wspace=0.4, hspace=0.4)
    plt.plot(t, R,label="расстояние от Солнца до Земли")
    plt.xlabel("t", fontsize=20)
    plt.ylabel("r", fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid(True)
    plt.legend(loc='lower right', bbox_to_anchor=(1, 0), 
           ncol=1, fontsize=12, framealpha=1)
    plt.savefig(r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\radius_vector.png")
    
    plt.figure(7)

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
            # Чтение времени
            line = stream.readline()
            if not line:
                break
            t = float(line.strip())
            time.append(t)
            
            # Чтение данных для N объектов
            for i in range(N):
                parts = stream.readline().strip().split()
                if len(parts) != 6:
                    raise ValueError(f"Ожидалось 6 значений для объекта {i}, получено {len(parts)}")
                
                x = float(parts[0])
                y = float(parts[1])
                # Остальные 4 значения игнорируем (temp в C++ коде)
                X[i].append(x)
                Y[i].append(y)
    
    # Построение графика
    plt.figure(figsize=(10, 8))
    plt.axis('equal')
    
    # Начальная точка (золотая)
    plt.scatter(X1, Y1, s=100, c='gold', linewidths=5)
    
    # Траектории
    for i in range(N):
        plt.plot(X[i], Y[i], linewidth=1)
    
    # Конечные точки (зеленые)
    for i in range(1, N):
        if X[i] and Y[i]:  # Проверка на пустые списки
            plt.scatter([X[i][-1]], [Y[i][-1]], s=100, c='green', linewidths=5)
    
    # Настройки графика
    plt.subplots_adjust(left=0.2, right=0.9, bottom=0.15, top=0.9, wspace=0.4, hspace=0.4)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel("x", fontsize=20)
    plt.ylabel("y", fontsize=20)
    plt.grid(True)
    
    # Сохранение графика
    output_path = r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\conversation_laws\traectories1.png"
    plt.savefig(output_path)
    plt.close()
    print(f"График сохранен: {output_path}")

def vs():
    with open(r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\positions.txt", 'r') as stream:
        N = int(stream.readline())
        time = []
        X = [0] * N
        Y = [0] * N
        j = 0
        
        while True:
            line = stream.readline()
            if not line:
                break
                
            parts = line.strip().split()
            if len(parts) != 1 + N*6:
                continue
                
            t = float(parts[0])
            time.append(t)
            
            for i in range(N):
                offset = 1 + i*6
                X[i] = float(parts[offset])
                Y[i] = float(parts[offset+1])
            
            plt.grid(True)
            plt.scatter(X, Y)
            plt.savefig(r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\positions\pos_eq" + str(j) + ".png")
            plt.clf()
            j += 1

def vectors_velocity():
    with open(r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\positions.txt", 'r') as stream:
        N = int(stream.readline())
        time = []
        X = [0] * N
        Y = [0] * N
        VX = [0] * N
        VY = [0] * N
        j = 0
        
        while True:
            line = stream.readline()
            if not line:
                break
                
            parts = line.strip().split()
            if len(parts) != 1 + N*6:
                continue
                
            t = float(parts[0])
            time.append(t)
            
            for i in range(N):
                offset = 1 + i*6
                X[i] = float(parts[offset])
                Y[i] = float(parts[offset+1])
                VX[i] = float(parts[offset+3])
                VY[i] = float(parts[offset+4])
            
            plt.grid(True)
            plt.quiver(X, Y, VX, VY)
            plt.savefig(r"C:\Users\mesho\Desktop\научка_2025_весна\программная_реализация_Равновесная_Модель\визуальзация_измерений\positions\vel_eq" + str(j) + ".png")
            plt.clf()
            j += 1