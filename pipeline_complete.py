import pygame
import numpy as np
import math
import copy
import time

pygame.init()

screen = pygame.display.set_mode((1280, 720))
clock = pygame.time.Clock()
running = True

def clamp(n, min, max):
    if n < min:
        return min
    elif n > max:
        return max
    else:
        return n

zfar = 600
znear = 0
fov = 0.78
tfov = math.tan(fov/2)
ra = 100/100
proj_matrix = np.matrix([ [1/((tfov)*ra), 0, 0, 0] , [0, 1/(tfov), 0, 0] , [0, 0, (zfar+znear)/(znear-zfar), 1] , [0, 0, (2*zfar*znear)/(znear - zfar), 0] ])
proj_matrix = proj_matrix.T

class Cuadrado:
    def __init__(self, pos_x, pos_y, tamano):
        self.rect = pygame.Rect(pos_x, pos_y, tamano, tamano)
    def draw(self, color=(64, 64, 64)):
        pygame.draw.rect(screen, color, self.rect, -1)
    def draw_fill(self, color="red"):
        pygame.draw.rect(screen, color, self.rect, 0)
class Grid:
    def __init__(self, filas, columnas, tamano_cuadradito, p0, p1):
        self.mapa = []
        self.color_buffer = []
        self.depth_buffer = []
        self.triangles_screen = []
        self.to_interpolate = []
        self.border = set()
        self.inside = set()
        self.fragments = {}
        self.transformed_fragments = {}
        self.p0 = p0
        self.p1 = p1
        self.filas = filas
        self.columnas = columnas
        self.reset_fragments()
        centr = pygame.Vector2(screen.get_width() / 2, screen.get_height() / 2)
        p0 = pygame.Vector2(centr.x - (tamano_cuadradito * columnas)/2, centr.y - (tamano_cuadradito * filas)/2)
        
        for f in range(0, filas):
            self.mapa.append([])
            self.color_buffer.append([])
            self.depth_buffer.append([])
            pos_y = p0.y + f*tamano_cuadradito
            
            for c in range(0, columnas):

                pos_x = p0.x + c*tamano_cuadradito
                
                self.mapa[-1].append(Cuadrado(pos_x, pos_y, tamano_cuadradito))
                self.color_buffer[-1].append((0,0,0))
                self.depth_buffer[-1].append(zfar)
    def calculate_todraw(self, p0, p1):
        # Algoritmo DDA

        # 1) Se calcula el n, Cant. de filas o columnas entre cada punto discreto.
        # Para contiguidad, nos aseguramos que haya si o si un punto por cada fila en caso de que haya más filas.
        # Caso contrario se eligen las columnas.

        n = max(abs(round(p0[0] - p1[0])), abs(round(p0[1] - p1[1]))) # Me devuelve o el máximo de filas, o el máximo de columnas entre c/u
        if (n == 0):
            return
        # Calculamos ambas derivadas, son el ratio de cambio continuo, el otro eje siempre se moverá en unidades de 1, 
        # mientras que este irá incrementando poco a poco según la "pendiente" o la dirección entre ambos puntos.
        dy = (p0[1] - p1[1])/n
        dx = (p0[0] - p1[0])/n
        (px, py, pz) = tuple(p1)

        # Simplemente nos fijamos cual gana.
        if abs(p0[0] - p1[0]) > abs(p0[1] - p1[1]):
            # Calculamos el incremento, si va pa un lado o pal otro
            incr = 1 if p0[0] - p1[0] > 0 else -1
            for i in range(0, n):
                # finalmente, arrancando de p0 le vamos sumando el incremento (1 o -1)
                px = round(px) + incr
                # al py le sumamos la pequeña derivada
                py = py + dy
                # agregamos el punto con un redondeo para que nos lo tire a n entero.
                if (round(py) > 0 and round(px) < self.columnas and round(py) < self.filas):
                    self.color_buffer[round(py)][round(px)] = (255, 255, 255)
        else:
            incr = 1 if p0[1] - p1[1] > 0 else -1
            for i in range(0, n):
                py = round(py) + incr
                px = px + dx

                if (round(py) > 0 and round(px) < self.columnas and round(py) < self.filas):
                    self.color_buffer[round(py)][round(px)] = (255, 255, 255)

    def send_triangles(self, triangles, extra_info):
        # triangles is a list of tuples with 3 elements, each being a 3D vertex.
        self.triangles = triangles
        self.extra_info = extra_info

    def procesamiento_vertices(self):
        self.transform()
    
    def project(self):
        cx = len(self.mapa[1])/2
        cy = len(self.mapa)/2

        self.triangles_screen = []
        self.to_interpolate = []
        eps = 1e-5
        for t in self.triangles:
            self.triangles_screen.append([])
            self.to_interpolate.append([])
            for p in t:

                pnp = np.array([p[0], p[1], p[2], 1])
                pnp.shape = [4, 1]
                pp = proj_matrix*pnp
                pp = pp * (1/pp[3, 0])
                pixel_x = cx*pp[0, 0] + cx
                pixel_y = cy * pp[1, 0] + cy
                self.triangles_screen[-1].append((pixel_x, pixel_y, pp[2, 0]))

    def calculate_to_draw_triangle(self):

        self.border = set()
        for t in self.triangles_screen:
            self.calculate_todraw(t[0], t[1])
            self.calculate_todraw(t[1], t[2])
            self.calculate_todraw(t[2], t[0])

    def reset_fragments(self):
        i = 0
        for f in self.mapa:
            j = 0
            for e in f:
                self.fragments[(i, j)] = []
                j+=1
            i+=1
    def reset_buffers(self):
        for i in range(0, self.filas):
            for j in range(0, self.columnas):
                self.color_buffer[i][j] = (0,0,0)
                self.depth_buffer[i][j] = zfar
    def interpolate_fragments_and_fill(self):
        triangle_counter = 0
        # Phong stuff
        luz_entrada = 0.75
        kambiente = np.array([50, 0, 50])
        kdifusa = np.array([255, 0, 255])
        kespecular = np.array([255, 255, 255])
        camera_dir = np.array([0, 0, 1])
        light_dir = np.array([1/math.sqrt(2), 1/math.sqrt(2), 0])
        q = 10
        for t in self.triangles_screen:
            p0 = np.array([ t[0][0] , t[0][1] , 0 ])
            p1 = np.array([ t[1][0] , t[1][1] , 0 ])
            p2 = np.array([ t[2][0] , t[2][1] , 0 ])

            p0_ = np.array([ t[0][0] , t[0][1] , self.triangles[triangle_counter][0][2] ])
            p1_ = np.array([ t[1][0] , t[1][1] , self.triangles[triangle_counter][1][2] ])
            p2_ = np.array([ t[2][0] , t[2][1] , self.triangles[triangle_counter][2][2] ])

            v01 = p0 - p1
            v02 = p2 - p0
            v21 = p1 - p2

            v01_ = p0_ - p1_
            v02_ = p2_ - p0_
            v21_ = p1_ - p2_

            at = np.cross(v01, v02)
            
            
            if (np.dot(at,at)==0):
                print("at ", at)
                print("p0 ", p0)
                print("p1 ", p1)
                print("p2", p2)
                print("v01 ", v01)
                print("v02 ", v02)
            at_ = np.cross(v01_, v02_)
            n = at_ / np.linalg.norm(at_)
            # Blinn - Phong
            h = (light_dir + camera_dir)/np.linalg.norm(light_dir + camera_dir)

            Id = kdifusa* ( np.dot(light_dir, n))
            Is = kespecular * ( np.dot(h, n)**q)

            Ir = luz_entrada * (kambiente + Id + Is)
            
            
            Ir = ( clamp(int(Ir[0]), 0, 255) , clamp(int(Ir[1]), 0, 255) , clamp(int(Ir[1]), 0, 255) )

            minx = max(min([p0[0], p1[0], p2[0]]), 0)
            maxx = min(max([p0[0], p1[0], p2[0]]), len(self.mapa[0])-1)

            miny = max(0, min([p0[1], p1[1], p2[1]]))
            maxy = min(max([p0[1], p1[1], p2[1]]), len(self.mapa)-1)

            

            for i in range(int(miny), int(maxy) + 1):
                for j in range(int(minx), int(maxx) + 1):
                    p = np.array([ j + 0.5 , i + 0.5 , 0 ])
                    a0 = np.cross(v21, p-p2)
                    a1 = np.cross(v02, p-p0)
                    a2 = np.cross(v01, p-p1)

                    alfa0 = np.dot(a0, at)/np.dot(at, at)
                    alfa1 = np.dot(a1, at)/np.dot(at, at)
                    alfa2 = np.dot(a2, at)/np.dot(at, at)

                    if (alfa0 < 0 or alfa1 < 0 or alfa2 < 0):
                        continue
                    
                    # hacemos un z-test del fragmento

                    # acordemonos que estos Z están siendo divididos por W.
                    zp0 = self.triangles[triangle_counter][0][2]
                    zp1 = self.triangles[triangle_counter][1][2]
                    zp2 = self.triangles[triangle_counter][2][2]

                    # 1) Interpolamos el inverso de Z: w = 1/z


                    # Por alguna razón z no se interpola linealmente pero el inverso si xd
                    wi = alfa0*(1/zp0) + alfa1*(1/zp1) + alfa2*(1/zp2)

                    zi = 1/wi
                    #print(zi)

                    if zi > self.depth_buffer[i][j]:
                        # Se descarta el fragmento
                        continue
                    else:
                        self.depth_buffer[i][j] = zi

                    # interpolación hiperbólica

                    # 1) Interpolamos 1/z

                    # 2) interpolamos el atributo peeero dividido por el Z de su punto correspondiente
                    r_color_p = alfa0*(255/zp0) + alfa1*(0/zp1) + alfa2*(0/zp2)
                    g_color_p = alfa1*(255/zp1)
                    b_color_p = alfa2*(255/zp2)

                    # Corregimos perspectiva
                    r_color = r_color_p / wi
                    g_color = g_color_p / wi
                    b_color = b_color_p / wi

                    # A VER: INTERPOLACIÓN HIPERBÓLICA: NO SE INTERPOLA EL ATRIBUTO, SE INTERPOLA EL ATRIBUTO DIVIDIDO POR Z
                    # PORQUE DE ESTA MANERA LA INTERPOLACION ES LINEAL.

                    # LUEGO, O ANTES, CUANDO SEA, INTERPOLAMOS LA INVERSA DE W EN ESE PUNTO, ENTONCES AL DIVIDIR POR W CANCELAMOS
                    # LA OPERACION PORQUE SE CANCELAN, NADA MAS QUE AHORA SI ESTA BIEN PORQUE EL W SI SE INTERPOLA BIEN,
                    # ES DECIR, LA INVERSA DE Z SI SE INTERPOLA LINEALMENTE. 

                    # ES UN TRUQUITO PARA CANCELAR LOS Z'S DE MANERA QUE QUEDE BIEN INTERPOLADO.

                    self.color_buffer[i][j] = Ir
            triangle_counter+=1
    def model_transform_to_origin(self):
        i = 0
        for t in self.triangles:
            j = 0
            for p in t:
                pass
                v = np.array([p[0], p[1], p[2], 1])

                v.shape = (4, 1)
                tm = np.matrix([[1, 0, 0, 0], (0, 1, 0, 0), (0, 0, 1, 0), (0, -80, -330, 1)])
                tm = tm.T

                rm = np.matrix([[math.cos(angulo), 0, -math.sin(angulo), 0], (0, 1, 0, 0), (math.sin(angulo), 0, math.cos(angulo), 0), (0, 0, 0, 1)])
                rm = rm.T

                btm = np.matrix([[1, 0, 0, 0], (0, 1, 0, 0), (0, 0, 1, 0), (0, 80, 330, 1)])
                btm = btm.T

                v = (tm)*v

                self.triangles[i][j] = (v[0, 0], v[1, 0], v[2, 0])
                j+=1
            i+=1
    def model_transform(self):
        i = 0
        for t in self.triangles:
            j = 0
            for p in t:
                pass
                v = np.array([p[0], p[1], p[2], 1])

                v.shape = (4, 1)
                tm = np.matrix([[1, 0, 0, 0], (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, -200, 1)])
                tm = tm.T

                rm = np.matrix([[math.cos(angulo), 0, -math.sin(angulo), 0], (0, 1, 0, 0), (math.sin(angulo), 0, math.cos(angulo), 0), (0, 0, 0, 1)])
                rm = rm.T

                btm = np.matrix([[1, 0, 0, 0], (0, 1, 0, 0), (0, 0, 1, 0), (0, 80, 330, 1)])
                btm = btm.T

                v = (btm)*v

                self.triangles[i][j] = (v[0, 0], v[1, 0], v[2, 0])
                j+=1
            i+=1
    def transform(self, angulo=0, altura=0):
        i = 0
        for t in self.triangles:
            j = 0
            for p in t:
                pass
                v = np.array([p[0], p[1], p[2], 1])

                v.shape = (4, 1)
                tm = np.matrix([[1, 0, 0, 0], (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, -200, 1)])
                tm = tm.T

                rm = np.matrix([[math.cos(angulo), 0, -math.sin(angulo), 0], (0, 1, 0, 0), (math.sin(angulo), 0, math.cos(angulo), 0), (0, 0, 0, 1)])
                rm = rm.T

                btm = np.matrix([[1, 0, 0, 0], (0, 1, 0, 0), (0, 0, 1, 0), (0, 80 + altura, 350, 1)])
                btm = btm.T

                v = (btm*rm)*v

                self.triangles[i][j] = (v[0, 0], v[1, 0], v[2, 0])
                j+=1
            i+=1
    def casteljau(self, puntos_control):
        divisions = 5

        delta_alfa = 1/divisions
        alfa = 0
        points = []
        for i in range(0, divisions + 1):
            print("Puntos de Control: ", puntos_control)
            casteljau_segments = self.build_segments_s(puntos_control)
            print("Segmentos: ", casteljau_segments)
            point = self.get_casteljau_point(casteljau_segments, alfa)
            points.append(point)

            alfa += delta_alfa
        return points

    def casteljau_surface(self, control_points):
        # Ahora control_points será una matriz de puntos en el espacio 
        # Ejemplo:
        # [ [ p0 , p1 , p2 , p3 ] 
        #   [ p4 , p5 , p6 , p7 ] ]

        points_matrix = []
        for i in range(0, len(control_points)):
            print(control_points[i])
            points_matrix.append(self.casteljau(control_points[i]))
            # points matrix tendrá la forma:
            # [ [i0 i1 i2 i3 i4 ]
            #   [i5 i6 i7 i8 i9 ]]
            # donde la primera fila son los puntos de la curva de bézier de la primer fila de puntos de control.  
            # y así.

        # Ahora aplicamos Casteljau de curvas de nuevo pero entre estos puntos, es decir todas las columnas.
        points_matrix_final = []
        for j in range(0, len(points_matrix[0])):
            # Quiero iterar las columnas y mandarle la columna entera a casteljau
            col = j  # índice de la columna (0-based)
            columna = [fila[col] for fila in points_matrix]
            points_matrix_final.append(self.casteljau(columna))

            # points matrix final tendría la forma
            # [ [f0 f1 f2 f3 f4]
            #   [f5 f6 f7 f8 f9]
            #   ...            ]

        # Quiero armar triángulos cruzados, tenemos que ir loopeando la matriz final.
        triangles = []
        for i in range(0, len(points_matrix_final) -1):
            for j in range(0, len(points_matrix_final[0]) - 1):
                t0 = [points_matrix_final[i][j], points_matrix_final[i][j+1], points_matrix_final[i+1][j+1]]
                t1 = [points_matrix_final[i][j], points_matrix_final[i+1][j+1], points_matrix_final[i+1][j]]
                triangles.append(t0)
                triangles.append(t1)
        return triangles

    def build_segments_s(self, points):
        segments = []
        for i in range(0, len(points) - 1):
            segments.append([])
            segments[-1].append(points[i])
            segments[-1].append(points[i+1])
        return segments
    def reset_points(self):
        self.points = []
    def get_casteljau_point(self, casteljau_segments, alfa):

        if (len(casteljau_segments) == 1):
            px = (1-alfa)*casteljau_segments[0][0][0] + alfa*casteljau_segments[0][1][0]
            py = (1-alfa)*casteljau_segments[0][0][1] + alfa*casteljau_segments[0][1][1]
            pz = (1-alfa)*casteljau_segments[0][0][2] + alfa*casteljau_segments[0][1][2]
            return (px, py, pz)

        pointss = []
        for i in range(0, len(casteljau_segments)):
            px = (1-alfa)*casteljau_segments[i][0][0] + alfa*casteljau_segments[i][1][0]
            py = (1-alfa)*casteljau_segments[i][0][1] + alfa*casteljau_segments[i][1][1]
            pz = (1-alfa)*casteljau_segments[0][0][2] + alfa*casteljau_segments[0][1][2]
            pointss.append( (px, py, pz) )
        casteljau_segments1 = self.build_segments_s(pointss)
        return self.get_casteljau_point(casteljau_segments1, alfa)
    def draw(self):
        for i in range(0, self.filas):
            for j in range(0, self.columnas):
                self.mapa[i][j].draw_fill(self.color_buffer[i][j])

grid = Grid(100, 100, 7, (5, 18), (19, 11))

extra_info = {"color": [(255, 0, 0) , (0, 255, 0) , (0, 0, 255)]}


control_points = [ [(100, 0, -100) , (50, -75, -100) , (-50, -75, -100) , (-100, 0, -100)],
                   [ (100, 0, 0) , (50, -30, 0) , (-50, -30, 0) , (-100, 0, 0) ],
                   [ (100, 0, 100) , (50, -75, 100) , (-50, -75, 100) , (-100, 0, 100) ]  ]

triangles = grid.casteljau_surface(control_points)
angulo = 0
altura = 0
wireframe = False

grid.send_triangles(copy.deepcopy(triangles), extra_info)

while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
    pressed_keys = pygame.key.get_pressed()
    modifiers_bitmask = pygame.key.get_mods()

    if pressed_keys[pygame.K_LEFT]:
        angulo+=0.1
        print(angulo)
    if pressed_keys[pygame.K_RIGHT]:
        angulo-=0.1
        print(angulo)
    if pressed_keys[pygame.K_UP]:
        altura-=3
    if pressed_keys[pygame.K_DOWN]:
        altura+=3

    screen.fill("black")
    grid.reset_buffers()

    grid.send_triangles(copy.deepcopy(triangles), extra_info)

    grid.transform(angulo, altura)
    
    if wireframe:
        grid.project()
        grid.calculate_to_draw_triangle()
    else:
        # Luego se rasteriza
        grid.project()
        # Interpolación
        grid.interpolate_fragments_and_fill()

    grid.draw()

    pygame.display.flip()

    dt = clock.tick(60) / 1000

pygame.quit()