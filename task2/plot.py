import matplotlib.pyplot as plt
import numpy as np

def plot_from_file(filename):
    """
    Рисует график функции по данным из файла с осью X через ноль оси Y
    """
    try:
        # Читаем данные из файла
        data = np.loadtxt(filename)
        x = data[:, 0]  # первая колонка - x
        y = data[:, 1]  # вторая колонка - y(x)
        
        # Создаем график
        plt.figure(figsize=(10, 6))
        plt.plot(x, y, 'b-', linewidth=2, label='S(x)')
        plt.grid(True, alpha=0.3)
        plt.xlabel('x', fontsize=12)
        plt.ylabel('y(x)', fontsize=12)
        plt.title('График функции', fontsize=14)
        plt.legend()
        
        # Добавляем ось X через ноль оси Y
        plt.axhline(y=0, color='black', linewidth=0.8)  # горизонтальная линия через y=0
        plt.axvline(x=0, color='black', linewidth=0.8)  # вертикальная линия через x=0
        
        # Показываем график
        plt.tight_layout()
        plt.show()
        
    except FileNotFoundError:
        print(f"Ошибка: Файл '{filename}' не найден")
    except Exception as e:
        print(f"Ошибка при чтении файла: {e}")

def plot_from_file_zero_axis(filename, title="График функции", xlabel="x", ylabel="y(x)"):
    """
    Рисует график функции с осью X через ноль оси Y и дополнительными настройками
    """
    try:
        # Читаем данные из файла
        data = np.loadtxt(filename)
        x = data[:, 0]
        y = data[:, 1]
        
        # Создаем график с настройками
        plt.figure(figsize=(12, 7))
        plt.plot(x, y, 'r-', linewidth=1.5, label=ylabel)
        
        # Настройки графика
        plt.grid(True, alpha=0.3)
        plt.xlabel(xlabel, fontsize=12)
        plt.ylabel(ylabel, fontsize=12)
        plt.title(title, fontsize=14)
        plt.legend(fontsize=10)
        
        # Добавляем оси через ноль
        plt.axhline(y=0, color='black', linewidth=1, linestyle='-')  # ось X
        plt.axvline(x=0, color='black', linewidth=1, linestyle='-')  # ось Y
        
        # Настраиваем отображение осей
        ax = plt.gca()
        
        # Перемещаем оси чтобы они проходили через ноль
        ax.spines['left'].set_position('zero')
        ax.spines['bottom'].set_position('zero')
        
        # Убираем правую и верхнюю оси
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        
        # Добавляем стрелки на концах осей
        ax.plot(1, 0, ">k", transform=ax.get_xaxis_transform(), clip_on=False)
        ax.plot(0, 1, "^k", transform=ax.get_yaxis_transform(), clip_on=False)
        
        # Добавляем сетку и улучшаем внешний вид
        plt.minorticks_on()
        plt.grid(which='minor', alpha=0.2)
        
        # Показываем график
        plt.tight_layout()
        plt.show()
        
        # Дополнительная информация о данных
        print(f"Количество точек: {len(x)}")
        print(f"Диапазон x: [{x.min():.3f}, {x.max():.3f}]")
        print(f"Диапазон y: [{y.min():.3f}, {y.max():.3f}]")
        
    except Exception as e:
        print(f"Ошибка: {e}")

def plot_classic_zero_axis(filename):
    """
    Классический вариант с осью X через ноль оси Y
    """
    try:
        data = np.loadtxt(filename)
        x = data[:, 0]
        y = data[:, 1]
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Рисуем график
        ax.plot(x, y, 'b-', linewidth=2)
        
        # Настраиваем оси чтобы проходили через ноль
        ax.spines['left'].set_position('zero')
        ax.spines['bottom'].set_position('zero')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        
        # Добавляем стрелки
        ax.plot(1, 0, ">k", transform=ax.get_xaxis_transform(), clip_on=False)
        ax.plot(0, 1, "^k", transform=ax.get_yaxis_transform(), clip_on=False)
        
        # Сетка
        ax.grid(True, alpha=0.3)
        ax.set_xlabel('x', fontsize=12)
        ax.set_ylabel('y(x)', fontsize=12)
        ax.set_title('График функции с осями через ноль', fontsize=14)
        
        plt.tight_layout()
        plt.show()
        
    except Exception as e:
        print(f"Ошибка: {e}")

def plot_simple_zero_axis(filename):
    """
    Простой вариант с осью X через ноль
    """
    data = np.loadtxt(filename)
    x = data[:, 0]
    y = data[:, 1]
    
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, 'b-', linewidth=1.5)
    
    # Добавляем ось X через y=0
    plt.axhline(y=0, color='k', linewidth=1)
    
    # Добавляем ось Y через x=0 если она в пределах графика
    if x.min() <= 0 <= x.max():
        plt.axvline(x=0, color='k', linewidth=1)
    
    plt.grid(True, alpha=0.3)
    plt.xlabel('x')
    plt.ylabel('y(x)')
    plt.title('График функции')
    plt.show()

# Пример использования
if __name__ == "__main__":
    # Простой вызов с осью X через ноль
    # plot_from_file("spline_output.txt")
    
    # Более продвинутый вариант с осями через ноль
    plot_from_file_zero_axis("out.txt", 
                           title="Кубический сплайн",
                           xlabel="Координата x", 
                           ylabel="S(x)")
    
    # Классический вариант
    # plot_classic_zero_axis("spline_output.txt")
    
    # Самый простой вариант
    # plot_simple_zero_axis("spline_output.txt")