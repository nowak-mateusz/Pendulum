#include <algorithm>
#include <cmath>
#include <map>
#include <iostream>
#include <vector>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <iostream>

// GLM Mathemtics
#include <glm/glm.hpp>

using namespace std;

//! Prezentacja i obliczenia dla whadeł sprzężonych
namespace pendulum
{
    const float PI = 3.14159265f; //!< stała PI
    const float DEG = PI / 180.f;  //!< 1 stopień [rad]

                                   //! zmodyfikowana funkcja sinus by poprawnie rysowały się elemty w układie współrzęnych OpenGL
    double fixsin(double rad)
    {
        return sin(PI + rad);
    }

    //! zmodyfikowana funkcja cosinus by poprawnie rysowały się elemty w układie współrzęnych OpenGL
    double fixcos(double rad)
    {
        return cos(PI + rad);
    }

    typedef glm::vec2 Point;             //!< Typ Punkt
    typedef glm::vec3 Color;             //!< Typ Kolor
    typedef std::vector<Point> Points;   //!< Typ kolekcji Punktów

    //! Wahadło
    /*!
    Reprezentacja graficzna pojedynczego wahadała.
    
    Położenie punktu masy wahadła obliczne jest za pomocą  pomocą równań:
    \f[x = px + l * sin(\alpha)\f]\f[y = py + l * cos(\alpha)\f]
    gdzie: \f$px\f$,\f$px\f$ - punkt zawieszenia wahadła,
    \f$l\f$ - długość nici
    \f$\alpha\f$ - kąt wahadła
    */

    class Pendulum
    {
    public:
        //! Konstruktor
        /*!
        \param p punkt "zawieszenia" wahadła w przestrzeni OpenGL
        \param length długość nici wahadła
        \param color kolor punktu masy
        */
        Pendulum(const Point& p, float length, const Color& color)
            : p0(p)
            , p1{ p.x, p.y + length }
            , l(length)
            , color(color)
        {}

        //! Konstruktor
        /*!
        \param x współrzędna x "zawieszenia" wahadła w przestrzeni OpenGL
        \param y wspołrzędna y "zawieszenia" wahadła w przestrzeni OpenGL
        \param length długość nici wahadła
        \param color kolor punktu masy
        */
        Pendulum(float x, float y, float length, const Color& color)
            : p0{ x, y }
            , p1{ x, y + length }
            , l(length)
            , color(color)
        {}

        //! Metoda obliczająca wspołrzęne potrzebne do narysowania wachadła
        /*!
        \param radian kąt wahadła w radianch
        */
        void calculate(double radian)
        {
            double sinVal = fixsin(radian);
            double cosVal = fixcos(radian);

            // wyznaczmy położenie punktu masy
            p1.x = p0.x + l * sinVal;
            p1.y = p0.y + l * cosVal;

            // wyznaczmy wspołrzędne graficznej reprezentacji punktu masy (bob)
            float cx = bobSize.x * cosVal;
            float sx = bobSize.x * sinVal;
            float cy = bobSize.y * cosVal;
            float sy = bobSize.y * sinVal;

            points.clear();
            points.push_back(Point{ p1.x - cx     , p1.y + sx });
            points.push_back(Point{ p1.x - cx + sy, p1.y + sx + cy });
            points.push_back(Point{ p1.x + cx + sy, p1.y - sx + cy });
            points.push_back(Point{ p1.x + cx     , p1.y - sx });
        }

        //! Metoda rysująca wahadło w OpenGL
        void draw()
        {
            drawLine(p0, p1);
            drawPolygon(points, color);
        }

        //! Funkcja pobierająca wspołrzędne graficznej reprezentacji punktu masy (bob)
        Points getBobPoints() const
        {
            return points;
        }

    private:

        //! Metoda rysująca linie prostą
        /*!
        \param p1 punkt początkowy linii
        \param p2 punkt końcowy linii
        */
        void drawLine(const Point& p1, const Point& p2)
        {
            glBegin(GL_LINES);
            glVertex2f(p1.x, p1.y);
            glVertex2f(p2.x, p2.y);
            glEnd();
        }

        //! Metoda rysująca wielokąt
        /*!
        \param point kolekcja punktów wielokąta
        \param color kolor wielokąta
        */
        void drawPolygon(const Points& points, const Color& color)
        {
            glColor3f(color.r, color.g, color.b);
            glBegin(GL_POLYGON);
            for (const Point& p : points)
                glVertex2f(p.x, p.y);
            glEnd();
            glColor3f(1, 1, 1);
        }

        const Point p0{ 0.f, 0.f };          //!< punkt zaczepienia wahadła
        Point p1{ p1.x, p1.y + l };          //!< punkt masy wahadła
        const float l{ 10.f };                //!< długośc rysowanego wahadła
        const Point bobSize{ 1.f, 3.f };      //!< rozmiar grafcznej reprezentacji punktu masy
        Points points;                        //!< punkty reprezentujące wilokąt dla graficznej reprezentacji punktu masy
        const Color color{ 1.f, 1.f, 1.f };   //!< kolor wahadła
    };

    //! Sprężyna
    /*!
    Reprezentacja graficzna sprężyny.
    \code{.unparsed}

    lHook  /\    /\    /\
    p0-----s0  \  /  \  /  \ s1-----p1
    \/    \/    \/

    \endcode
    */

    class Spring
    {
    public:

        //! Metoda obliczająca wspołrzęne potrzebne do narysowania sprężyny
        /*!
        \param x0 współrzędna x początku zaczepienia sprężyny
        \param y0 współrzędna y początku zaczepienia sprężyny
        \param x1 współrzędna x końca zaczepienia sprężyny
        \param y1 współrzędna y końca zaczepienia sprężyny
        */
        void calculate(float x0, float y0, float x1, float y1)
        {
            isDraw = true;

            // uaktualnienie punktów początkowych i końcowych sprężyny
            p0.x = x0;
            p0.y = y0;

            p1.x = x1;
            p1.y = y1;

            float dx = p1.x - p0.x;
            float dy = p1.y - p1.y;
            float l = sqrt(dx * dx + dy * dy); // całkowita długość sprężyny

            if (l < 2 * lHook)
            {
                isDraw = false;
                return;
            }

            float q = lHook / l; // początek spirali jako ułamek całkowitej długości

            // uaktualnienie punktów początku i końca spirali
            s0.x = p0.x + q * dx;
            s0.y = p0.y + q * dy;
            s1.x = p1.x - q * dx;
            s1.y = p1.y - q * dy;

            float dsx = s1.x - s0.x;
            float dsy = s1.y - s0.y;
            l = sqrt(dsx * dsx + dsy * dsy); // długość spirali

            int maxNumCoilsPoints = numCoils * 360 / stepCoils;  // maksymalna ilość punktów spiralii

            coilsPoints.clear();

            // wyznacznie kolejnych punktów spirali
            for (int i = 1; i <= maxNumCoilsPoints; i++)
            {
                float a = i / (float)maxNumCoilsPoints;
                float b = (springHeight / l) * fixsin(i * stepCoils * DEG);

                Point coil;
                coil.x = s0.x + a * dsx + b * dsy;
                coil.y = s0.y + a * dsy - b * dsx;

                coilsPoints.push_back(coil);
            }
        }

        //! Metoda rysująca sprężynę w OpenGL
        void draw()
        {
            if (!isDraw)
                return;

            glBegin(GL_LINE_STRIP);
            glVertex2f(p0.x, p0.y);
            glVertex2f(s0.x, s0.y);

            for (const Point& p : coilsPoints)
                glVertex2f(p.x, p.y);

            glVertex2f(p1.x, p1.y);
            glEnd();
        }

    private:
        bool  isDraw{ true };             //!< czy możemy narysować sprężyne
        Point p0;                         //!< punkt początku sprężyny
        Point p1;                         //!< punkt końca sprężyny
        Point s0;                         //!< punkt początku spirali 
        Point s1;                         //!< punkt zakończenia spirali
        float lHook{ 1.f };               //!< długośc zaczepienia sprężyny
        const int   numCoils{ 10 };       //!< liczba zwojów
        const int   stepCoils{ 90 };      //!< wielkość kroku spirali (stopnie)
        const float springHeight{ 1.f };  //!< wysokość sprężyny 
        Points coilsPoints;               //!< kolekcja kolejnych punktów spiralii
    };

    //! Diagram
    /*!
    Klasa reprezentująca wykres
    */
    class Diagram
    {
    public:
        //! Konstruktor
        /*!
        \param x współrzędna x wstawienia Diagramu w przestrzeni OpenGL
        \param y wspołrzędna y wstawienia Diagramu w przestrzeni OpenGL
        \param color kolor wskażnika
        */
        Diagram(float x, float y, const Color& color, float omega1, float omega2)
            : p0{ x, y }
            , color(color)
            , omega1(omega1)
            , omega2(omega2)
        {}

        void setangle(float angle1, float angle2)
        {
            a1 = angle1;
            a2 = angle2;
        }

        //! Metoda obliczająca przebieg wykresu
        /*!
        \param t czas symulacji
        \param y0 współrzędna y początku zaczepienia sprężyny
        \param x1 współrzędna x końca zaczepienia sprężyny
        \param y1 współrzędna y końca zaczepienia sprężyny
        */
        void calculate(float t)
        {
            float pixT = size.x / (2 * dt); // współczynnik konwersji dla osi poziomej (pikseli na sekundę)
            float pixY = size.y;            // współczynnik konwersji dla osi pionowej
            float t0 = t - min(t, dt);      // cas początku wykresu
            float yy = a1 * fixcos(omega1*t0) + a2 * fixcos(omega2*t0); // odchylenie dla punktu początkowego
            float x = p0.x;                 // wspórzędna x dla punktu począkowego
            float y = p0.y - pixY * yy;     //wspołrzęna y dla punktu począkowego

                                            //uaktualnienie punktu początkowego
            d0.x = x;
            d0.y = y;

            points.clear();

            //oblicznie kolejnych punktów
            while (x < p0.x + size.x)
            {
                x += 0.01;
                float tt = (x - p0.x) / pixT + t0;
                yy = a1 * fixcos(omega1*tt) + a2 * fixcos(omega2*tt); //odchylania dla kolejnych próbek czasu
                y = p0.y - pixY * yy;

                Point p{ x, y };
                points.push_back(p);
            }

            yy = a1 * fixcos(omega1*t) + a2 * fixcos(omega2*t); //odchylania dla aktualnej chwili czasu

                                                                //wyznaczamy położenie wskażnika
            pointer.x = p0.x + min(size.x / 2, t*pixT);
            pointer.y = p0.y - pixY * yy;
        }

        void draw()
        {
            //rysujemy wykres
            glBegin(GL_LINE_STRIP);
            glVertex2f(d0.x, d0.y);

            for (const Point& p : points)
                glVertex2f(p.x, p.y);

            glEnd();

            //osie
            
            glBegin(GL_LINES);
            glVertex2f(p0.x, p0.y );
            glVertex2f(p0.x + size.x, p0.y);
            glEnd();

            glBegin(GL_LINES);
            glVertex2f(p0.x, p0.y- size.y/4);
            glVertex2f(p0.x, p0.y + size.y/4);
            glEnd();

            //rysujemy wskażnik
            drawCircle(pointer.x, pointer.y, 0.25, color);
        }

    private:
        //! metoda rysująca koło jako wielokąt
        /*!
        \param sx współrzędna x środka koła
        \param sy wspołrzędna y środka koła
        \param r promień koła
        \param color kolor koła
        */
        void drawCircle(float sx, float sy, float r, const Color& color)
        {
            const float segments = 25;
            float theta = 2 * PI / segments;
            float tan_factor = tanf(theta);
            float radial_factor = cosf(theta);
            float x = r;
            float y = 0;

            int cache_pt = 0;
            double cache_x;
            double cache_y;

            glColor3f(color.r, color.g, color.b);
            glBegin(GL_POLYGON);
            for (int ii = 0; ii < segments; ii++) {
                if (!cache_pt) {
                    cache_x = x + sx;
                    cache_y = y + sy;
                    cache_pt = 1;
                }
                else {
                    //glVertex2f(cache_x,cache_y);
                    glVertex2f(x + sx, y + sy);
                    cache_x = x + sx;
                    cache_y = y + sy;
                }
                float tx = -y;
                float ty = x;
                x += tx * tan_factor;
                y += ty * tan_factor;
                x *= radial_factor;
                y *= radial_factor;
            }
            glEnd();
            glColor3f(1.f, 1.f, 1.f);
        }

        const Point p0{ 0.f, 0.f };    //!< punkt wstawienia diagramu 
        Point d0{ 0.f, 0.f };          //!< punkt początku wykresu
        Point size{ 10.f, 10.f };      //!< szerokośći, wysokość wykresu
        Points points;                 //!< kolejne punkty na wykresie
        Point pointer{ 0, 0 };         //!< połorzenie wskażnika na wykresie
        const float dt{ 20.f };        //!< buffor czasu  
        const Color color;             //!< kolor znacznika
        const float omega1;
        const float omega2;
        float a1;
        float a2;
    };


    //! Symulacja
    /*!
    Klasa reprezentująca symulacje. Wykonuje obliczenia i wyświetla opowiednie elementy.
    */
    class Simulation
    {
    public:
        //! Konstruktor
        Simulation()
            : pendulum1(-10.f, 10.f, 20.f, Color{ 1.f, 0.f, 0.f })
            , pendulum2(10.f, 10.f, 20.f, Color{ 0.f, 0.f, 1.f })
            , diagram1(15.f, 5.f, Color{ 1.f, 0.f, 0.f }, omega1, omega2)
            , diagram2(15.f, -2.f, Color{ 0.f, 0.f, 1.f }, omega1, omega2)
        {
            diagram1.setangle(a1, a2);
            diagram2.setangle(a1, -a2);
        }

        //! metoda pozwala ustawić kąty początkowe wachadeł
        /*!
        \param angle1 kąt wychylenia początkowego wahadła 1 [stopnie]
        \param angle2 kąt wychylenia początkowego wahadła 2 [stopnie]
        */
        void setStartAngle(float angle1, float angle2)
        {
            alpha01 = angle1 * DEG;
            alpha02 = angle2 * DEG;
            a1 = (alpha01 + alpha02) / 2.f;
            a2 = (alpha01 - alpha02) / 2.f;

            diagram1.setangle(a1, a2);
            diagram2.setangle(a1, -a2);
        }

        //! metoda wykonująca obliczenia dla podanej chwilii czasu
        /*!
        \param t czas symulacji dla którego mają być wykonane obliczenia
        */
        void calculate(float t)
        {
            float cos1 = fixcos(omega1*t);
            float cos2 = fixcos(omega2*t);
            float alpha1 = a1 * cos1 + a2 * cos2;
            float alpha2 = a1 * cos1 - a2 * cos2;

            pendulum1.calculate(alpha1);
            pendulum2.calculate(alpha2);

            Points polygon1 = pendulum1.getBobPoints();
            Points polygon2 = pendulum2.getBobPoints();

            float fx1 = (polygon1[0].x + polygon1[1].x) / 2;
            float fy1 = (polygon1[0].y + polygon1[1].y) / 2;
            float fx2 = (polygon2[2].x + polygon2[3].x) / 2;
            float fy2 = (polygon2[2].y + polygon2[3].y) / 2;
            spring.calculate(fx1, fy1, fx2, fy2);

            diagram1.calculate(t);
            diagram2.calculate(t);
        }

        //! metoda rysująca elemty symulacjii
        void draw()
        {
            pendulum1.draw();
            pendulum2.draw();
            spring.draw();
            diagram1.draw();
            diagram2.draw();
        }

    private:

        const float g{ 9.81f };  //!< \f$g\f$ stała przyspieszenia ziemskiego [m/s] 
        const float l{ 1.f };    //!< \f$l\f$ długość wahadła [m]
        const float k{ 1.f };    //!< \f$k\f$ stała sprężystości [N/m]
        const float m{ 1.f };    //!< \f$m\f$ masa ciała wahadła [kg]
                                 //float time { 0.f };       //!< czas symulacji
        const float omega1{ sqrt(g / l) };              //!< \f$\omega_1\f$  częstotliwość własna układu
        const float omega2{ sqrt(g / l + 2 * k / m) };  //!< \f$\omega_2\f$  częstotliwość własna układu
        float alpha01{ -10.f * DEG };                     //!< \f$\alpha_1\f$  kąt początkowy pierwszego wahadła [rad]
        float alpha02{ 0.f * DEG };                       //!< \f$\alpha_2\f$  kąt początkowy drugieg wahadła  [rad]
        float a1{ (alpha01 + alpha02) / 2.f };          //!< zmienna pomocnicza (max aplituda?)
        float a2{ (alpha01 - alpha02) / 2.f };          //!< zmienna pomocnicza (min aplituda?)

        Pendulum pendulum1;     //!< Wahadło 1
        Pendulum pendulum2;     //!< Wahadło 2
        Spring   spring;        //!< Sprężyna
        Diagram  diagram1;      //!< Wykres dla wahadła 1
        Diagram  diagram2;      //!< Wykres dla wahadła 2
    };
}

static  pendulum::Simulation simulation;  //!< globalny obiekt symulacjii
    
//// OpenGL stuff ////

//! metoda wysietlająca w OpenGL
void display(float time)
{
    glLoadIdentity();
    glTranslatef(-5, 5, -30);
    glColor3f(1, 1, 1);

    simulation.calculate(time);
    simulation.draw();
}

//! metoda wyznaczająca maczierz perspektywiczną 
void perspectiveGL(GLdouble fovY, GLdouble aspect, GLdouble zNear, GLdouble zFar)
{
    const GLdouble pi = 3.1415926535897932384626433832795;
    GLdouble fW, fH;

    //fH = tan( (fovY / 2) / 180 * pi ) * zNear;
    fH = tan(fovY / 360 * pi) * zNear;
    fW = fH * aspect;

    glFrustum(-fW, fW, -fH, fH, zNear, zFar);
}

//! metoda przeskalująca wywołana po zmainie szerokośći okna
void reshape(int w, int h)
{
    glViewport(0, 0, (GLsizei)w, (GLsizei)h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    perspectiveGL(60.0, (GLfloat)w / (GLfloat)h, 0.1, 100.0);
    glMatrixMode(GL_MODELVIEW);

    glLoadIdentity();
}

static int winWidth = 800;      //!< szerokość okna
static int winHeight = 600;     //!< wysokość okna
static bool  start  = false;    //!< czy symulacja jest uruchomina
static float time = 0.f;        //!< czas symulacji
static float timeoffset = 0.02f; //!< offset czasu symulacji
static float angle1 = 10;       //!< kąt początkowy pierwszego wahadła
static float angle2 = 0;        //!< kąt początkowy drugiego wahadła

//! metoda odpowiadająca za wyświetlanie przycisków
void displayGui()
{
    const int minWidth = 350;
    const int minHeight = 60;

    ImGui::SetNextWindowPos(ImVec2((winWidth - minWidth) / 2, (winHeight- minHeight - 20)));
    ImGui::SetNextWindowSize(ImVec2(minWidth, minHeight));
    ImGui::Begin("Ustawienia", NULL, ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove);

    ImGui::BeginGroup();
    ImGui::BeginGroup();
    if (ImGui::Button("Start", ImVec2(80.f, 40.f)))
        start = true;
    ImGui::SameLine();
    if (ImGui::Button("Stop", ImVec2(80.f, 40.f)))
        start = false;

    ImGui::EndGroup();
    ImGui::SameLine();
    ImGui::Spacing();
    ImGui::SameLine();
    ImGui::BeginGroup();
    ImGui::PushItemWidth(60.f);
    ImGui::Text("alpha1");
    ImGui::SameLine();
    float currentAngle1 = angle1;
    float currentAngle2 = angle2;
    if (ImGui::InputFloat("##alpha1", &angle1, 0.f, 0.f, 2, ImGuiInputTextFlags_EnterReturnsTrue))
    {
        if (start)
            angle1 = currentAngle1;

        if (angle1 < -10.f)
            angle1 = -10.f;
        if (angle1 > 10.f)
            angle1 = 10.f;

        if (!start)
        {
            simulation.setStartAngle(angle1, angle2);
            time = 0.f;
        }
    }
    ImGui::Text("alpha2");
    ImGui::SameLine();
    if (ImGui::InputFloat("##alpha2", &angle2, 0.f, 0.f, 2, ImGuiInputTextFlags_EnterReturnsTrue))
    {
        if (start)
            angle2 = currentAngle2;

        if (angle2 < -10.f)
            angle2 = -10.f;
        if (angle2 > 10.f)
            angle2 = 10.f;

        if (!start)
            simulation.setStartAngle(angle1, angle2);
    }
    ImGui::EndGroup();
    ImGui::EndGroup();
    ImGui::End();
}

int main(void)
{
    GLFWwindow* window;

    /// inicjalizacja biblioteki glfw
    if (!glfwInit())
        return -1;

    /// tworzenie okna i kontekstu dla OpenGL
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
    window = glfwCreateWindow(winWidth, winHeight, "Pendulum", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1); // Enable vsync

    GLenum err = glewInit();
    if (GLEW_OK != err)
        fprintf(stderr, "Error: %s\n", glewGetErrorString(err));

    fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));

    reshape(winWidth, winHeight);

    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); 
    (void)io;
    io.IniFilename = NULL;
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");
    ImGui::StyleColorsDark();

    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.f);
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.f,0.f,0.f,0.f));
    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window))
    {
        /* Render here */
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        display(time);
        displayGui();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events */
        glfwPollEvents();

        if (start)
            time += timeoffset;
    }

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}