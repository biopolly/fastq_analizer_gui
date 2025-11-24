import gzip
import os
import tkinter as tk
from tkinter import filedialog, messagebox
from collections import defaultdict




class FASTQProcessor:
    """класс для обработки fastq файлов"""
    
    def __init__(self):
        
        
        self.sequences_count = 0          # количество последовательностей
        self.total_length = 0             # общая длина всех последовательностей
        self.sequence_lengths = []        # список длин каждой последовательности
        self.quality_scores = []          # список оценок качества
        self.base_counts = defaultdict(int) # счетчик оснований
        
    def parse_fastq(self, file_path):
        """
        парсинг fastq файла и сбор статистики
        
        параметры:
            file_path путь к файлу для анализа
        """
        # сброс предыдущей статистики
        self.sequences_count = 0
        self.total_length = 0
        self.sequence_lengths = []
        self.quality_scores = []
        self.base_counts = defaultdict(int)
        
        # выбор способа открытия файла (обычный или gzip)
        if file_path.endswith('.gz'):
            # открытие gzip файла в текстовом режиме
            opener = gzip.open
        else:
            # открытие обычного текстового файла
            opener = open
        
        # открытие файла с обработкой ошибок кодировки
        with opener(file_path, 'rt', encoding='utf-8', errors='ignore') as file:
            # чтение всех строк файла
            lines = file.readlines()
            
        # обработка файла по блокам по 4 строки 
        
        for i in range(0, len(lines), 4):
            # проверка что есть достаточно строк для полного блока
            if i + 3 < len(lines):
                
                # строка с последовательностью 
                
                sequence_line = lines[i + 1].strip()
                quality_line = lines[i + 3].strip()
                
                # увеличение счетчика последовательностей
                self.sequences_count += 1
                # вычисление длины текущей последовательности
                seq_len = len(sequence_line)
                # добавление к общей длине
                self.total_length += seq_len
                # сохранение длины последовательности
                self.sequence_lengths.append(seq_len)
                
                # подсчет оснований в последовательности
                self._count_bases(sequence_line)
                
                self._analyze_quality(quality_line)

    
    def _count_bases(self, sequence):
        """
        подсчет количества каждого типа оснований в последовательности
        
        параметры sequence строка с последовательностью нуклеотидов
        """
       

        sequence_upper = sequence.upper()
        
        # подсчет каждого основания в последовательности
        for base in sequence_upper:
            if base in 'ATGC':
                self.base_counts[base] += 1
    
    def _analyze_quality(self, quality_string):
        """
        анализ строки с оценками качества
        
        параметры:
            quality_string строка с символами качества
        """
        # преобразование символов качества в числовые значения
        # fastq  кодируется как ord 33
        quality_scores = [ord(char) - 33 for char in quality_string]
        # добавление оценок качества в общий список
        self.quality_scores.extend(quality_scores)
    
    def get_basic_stats(self):
        """получение базовой статистики по файлу"""
        # проверка наличия данных
        if self.sequences_count == 0:
            return "нет данных для анализа"
            
        # вычисление средней длины последовательности
        avg_length = self.total_length / self.sequences_count
        
        # формирование строки со статистикой
        stats = f"общее количество последовательностей: {self.sequences_count}\n"
        stats += f"средняя длина последовательности: {avg_length:.1f} bp\n"
        stats += f"минимальная длина: {min(self.sequence_lengths)} bp\n"
        stats += f"максимальная длина: {max(self.sequence_lengths)} bp\n"
        stats += f"общая длина всех последовательностей: {self.total_length} bp\n"
        
        return stats
    
    def get_quality_stats(self):
        """получение статистики по качеству последовательностей"""
        if not self.quality_scores:
            return "нет данных о качестве"
            
        # вычисление статистик качества
        avg_quality = sum(self.quality_scores) / len(self.quality_scores)
        min_quality = min(self.quality_scores)
        max_quality = max(self.quality_scores)
        
        stats = f"среднее качество: {avg_quality:.2f}\n"
        stats += f"минимальное качество: {min_quality}\n"
        stats += f"максимальное качество: {max_quality}\n"
        stats += f"общее количество оценок качества: {len(self.quality_scores)}"
        
        return stats
    
    def get_base_composition(self):
        """получение состава оснований"""
        if not self.base_counts:
            return "нет данных о составе оснований"
            
        # вычисление общего количества оснований
        total_bases = sum(self.base_counts.values())
        
        # формирование строки с составом оснований
        composition = "состав оснований:\n"
        for base in 'ATGC':
            count = self.base_counts[base]
            percentage = (count / total_bases) * 100 if total_bases > 0 else 0
            composition += f"  {base}: {count} ({percentage:.1f}%)\n"
        
        return composition
    




class FASTQAnalyzer:
    """класс графического интерфейса для анализа fastq файлов"""
    
    def __init__(self, root):
        """
        инициализация интерфейса
        
        параметры:
            root - корневое окно приложения
        """
        self.root = root
        # установка заголовка окна
        self.root.title("анализатор fastq файлов")
        # установка начального размера окна
        self.root.geometry("700x500")
        # запрет изменения размера окна
        self.root.resizable(True, True)
        
    

        self.processor = FASTQProcessor()
        # переменная для хранения пути к текущему файлу
        self.current_file = None
        
        
        self.create_interface()
    
    def create_interface(self):
        """создание пользовательского интерфейса"""
        
        main_frame = tk.Frame(self.root, padx=20, pady=20)
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # заголовок приложения
        title_label = tk.Label(
            main_frame, 
            text="анализатор fastq файлов", 
            font=("Times New Roman", 20, "bold"),
            fg="darkgreen"
        )
        title_label.pack(pady=10)
        
        # фрейм для кнопок управления
        button_frame = tk.Frame(main_frame)
        button_frame.pack(pady=10)
        
        # кнопка выбора файла
        self.select_btn = tk.Button(
            button_frame, 
            text="выбрать fastq файл", 
            command=self.select_file, 
            width=20, 
            height=2,
            bg="lightyellow",
            font=("Arial", 10)
        )
        self.select_btn.pack(side=tk.LEFT, padx=5)
        
        # кнопка анализа файла 
        self.analyze_btn = tk.Button(
            button_frame, 
            text="проанализировать", 
            command=self.analyze_file, 
            state=tk.DISABLED,
            width=15, 
            height=2,
            bg="lightpink",
            font=("Arial", 10)
        )
        self.analyze_btn.pack(side=tk.LEFT, padx=5)
        
        # метка для отображения информации о выбранном файле
        self.file_label = tk.Label(
            main_frame, 
            text="файл не выбран", 
            wraplength=650, 
            justify="left",
            bg="white",
            relief="sunken",
            padx=10,
            pady=5
        )
        self.file_label.pack(pady=5, fill=tk.X)
        
        # создание текстовых областей для разных типов статистики
        self.create_text_areas(main_frame)
    
    def create_text_areas(self, parent):
        """создание текстовых областей для отображения статистики"""
        # фрейм для текстовых областей
        text_frame = tk.Frame(parent)
        text_frame.pack(pady=10, fill=tk.BOTH, expand=True)
        
        # создание текстовых виджетов с метками
        self.create_text_widget(text_frame, "базовая статистика", 0)
        self.create_text_widget(text_frame, "качество последовательностей", 1)
        self.create_text_widget(text_frame, "состав оснований", 2)
    
    def create_text_widget(self, parent, title, row):
        """
        создание текстового виджета с заголовком
        
        
        parent  родительский фрейм
        title заголовок текстовой области
        row номер строки для размещения
        """
        # создание фрейма для текстового виджета
        frame = tk.Frame(parent)
        frame.grid(row=row, column=0, sticky="ew", pady=5)
        frame.columnconfigure(0, weight=1)
        
        # метка с заголовком
        label = tk.Label(frame, text=title, font=("Arial", 11, "bold"))
        label.grid(row=0, column=0, sticky="w", pady=(0, 5))
        
        # текстовое поле для отображения статистики
        text_widget = tk.Text(
            frame, 
            height=6, 
            width=80, 
            font=("Courier", 9),
            wrap=tk.WORD
        )
        text_widget.grid(row=1, column=0, sticky="ew", padx=5)
        
        # добавление скроллбара для текстового поля
        scrollbar = tk.Scrollbar(frame, command=text_widget.yview)
        scrollbar.grid(row=1, column=1, sticky="ns")
        text_widget.config(yscrollcommand=scrollbar.set)
        
        # сохранение ссылки на текстовый виджет
        if row == 0:
            self.basic_text = text_widget
        elif row == 1:
            self.quality_text = text_widget
        elif row == 2:
            self.composition_text = text_widget
    
    def select_file(self):
        """выбор файла через диалоговое окно"""
        # определение фильтров для типов файлов
        file_types = [
            ("fastq files", "*.fastq *.fq"),
            ("gzipped fastq", "*.fastq.gz *.fq.gz"),
            ("all files", "*.*")
        ]
        
        # открытие диалогового окна выбора файла
        filename = filedialog.askopenfilename(
            title="выберите fastq файл",
            filetypes=file_types
        )
        
        # обработка выбранного файла
        if filename:
            self.current_file = filename
            # обновление метки с информацией о файле
            file_info = f"выбран файл: {os.path.basename(filename)}"
            self.file_label.config(text=file_info)
            # активация кнопки анализа
            self.analyze_btn.config(state=tk.NORMAL)
    
    def analyze_file(self):
        """анализ выбранного файла"""
        # проверка что файл выбран
        if not self.current_file:
            messagebox.showwarning("внимание", "сначала выберите файл")
            return
            
        try:
            # отображение сообщения о начале обработки
            self.show_loading_message()
            
            # анализ файла с помощью процессора
            self.processor.parse_fastq(self.current_file)
            
            # отображение результатов анализа
            self.display_results()
            
        except Exception as e:
            # обработка ошибок при анализе
            error_message = f"ошибка при обработке файла:\n{str(e)}"
            messagebox.showerror("ошибка", error_message)
            self.clear_results()
    
    def show_loading_message(self):
        """отображение сообщения о загрузке"""
        loading_text = "обработка файла...\nпожалуйста, подождите."
        
        # очистка всех текстовых полей и отображение сообщения о загрузке
        for text_widget in [self.basic_text, self.quality_text, self.composition_text]:
            text_widget.delete(1.0, tk.END)
            text_widget.insert(tk.END, loading_text)
        
        # обновление интерфейса
        self.root.update()
    
    def display_results(self):
        """отображение результатов анализа"""
        # получение статистики из процессора
        basic_stats = self.processor.get_basic_stats()
        quality_stats = self.processor.get_quality_stats()
        composition_stats = self.processor.get_base_composition()
        
        # получение информации о размере файла
        file_size = os.path.getsize(self.current_file)
        file_info = f"\n\nинформация о файле:\n"
        file_info += f"размер файла: {file_size / 1024 / 1024:.2f} MB\n"
        file_info += f"имя файла: {os.path.basename(self.current_file)}"
        
        # отображение базовой статистики
        self.basic_text.delete(1.0, tk.END)
        self.basic_text.insert(tk.END, basic_stats + file_info)
        
        # отображение статистики качества
        self.quality_text.delete(1.0, tk.END)
        self.quality_text.insert(tk.END, quality_stats)
        
        # отображение состава оснований
        self.composition_text.delete(1.0, tk.END)
        self.composition_text.insert(tk.END, composition_stats)
    
    def clear_results(self):
        
        # очистка всех текстовых полей
        self.basic_text.delete(1.0, tk.END)
        self.quality_text.delete(1.0, tk.END)
        self.composition_text.delete(1.0, tk.END)



def main():

    root = tk.Tk()
    # создание экземпляра приложения
    app = FASTQAnalyzer(root)
   
    root.mainloop()



if __name__ == "__main__":
    main()