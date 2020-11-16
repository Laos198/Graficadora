using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace XALIB
{
    public class Funciones
    { /// <summary>
    /// Suma de dos valores tipo double, devuelve el valor de la suma 
    /// </summary>
    /// <param name="A"></param>
    /// <param name="B"></param>
    /// <returns></returns>
        public double sum(double A, double B)
        {
            double r = A + B;
            return r;
        }
        /// <summary>
        /// Retorna la suma de una lista
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public double sum(List<double> x)
        {
            double suma = 0;
            for (int i = 0; i < x.Count(); i++)
            {
                suma = suma + x[i];
            }
            return suma;
        }
        /// <summary>
        /// Retorna el valor promedio de una lista de datos 
        /// </summary>
        /// <param name="datos"></param>
        /// <returns></returns>
        public double mean(List<double> datos)
        {
            double prom;
            double suma = 0;
            double n = datos.Count();
            for (int i = 0; i <= n - 1; i++)
            {
                suma = suma + datos[i];

            }
            prom = suma / n;
            return prom;
        }
        /// <summary>
        ///  retorna el valor mínimo de una lista de datos (Checar)
        /// </summary>
        /// <param name="datos"></param>
        /// <returns></returns>
        public double max(List<double> datos)
        {
            double maximo;
            List<double> datos2 = new List<double>();
            datos2 = datos;
            datos2.Sort();
            maximo = datos2[0];
            return maximo;
        }
        /// <summary>
        /// Retorna el calor máximo de una lista de datos 
        /// </summary>
        /// <param name="datos"></param>
        /// <returns></returns>
        public double min(List<double> datos)
        {
            double min;
            List<double> datos2 = new List<double>();
            datos2 = datos;
            datos2.Sort();
            min = datos2[datos2.Count - 1];
            return min;
        }
        /// <summary>
        ///  Retorna la varianza de una lista de datos 
        /// </summary>
        /// <param name="datos"></param>
        /// <returns></returns>
        public double var(List<double> datos)
        {
            double variance;
            double promedio = mean(datos);
            double n = datos.Count();
            double suma = 0;
            for (int i = 0; i <= n - 1; i++)
            {
                suma = suma + Math.Pow(datos[i] - promedio, 2);
            }
            variance = suma / (n);

            return variance;
        }
        /// <summary>
        /// std devuelve la desviacion estandar de un arreglo (Lista)
        /// </summary>
        /// <param name="datos"></param>
        /// <returns></returns>
        public double std(List<double> datos)
        {
            double desv_std = Math.Sqrt(var(datos));
            return desv_std;
        }
        /// <summary>
        ///  Genera un vector de x1 a x2 con numero de puntos n 
        /// </summary>
        /// <param name="x1"></param>
        /// <param name="x2"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public double[] linspace_2(double x1, double x2, int n)
        {
            double[] lin = new double[n];
            double paso = (x2 - x1) / (n - 1);
            for (int i = 0; i < n; i++)
            {
                lin[i] = x1 + paso * i;
            }
            return lin;
        }
        /// <summary>
        /// Te devuelve un vector de x1-x2 con n pasos 
        /// </summary>
        /// <param name="x1"></param>
        /// <param name="x2"></param>
        /// <param name="n"></param>
        /// <returns></returns>
        public double[] linspace(double x1, double x2, double n)
        {
            double paso = n;
            int l = Convert.ToInt32(1 + ((x2 - x1) / paso));

            double[] lin = new double[l];
            for (int i = 0; i < l; i++)
            {
                lin[i] = x1 + paso * i;
            }
            return lin;
        }

        /// <summary>
        /// Regresa un arreglo con la distribucion normal con parametros, una lista, el promedio (mu) y la desviación estandar (sigma)
        /// </summary>
        /// <param name="datos"></param>
        /// <param name="mu"></param>
        /// <param name="sigma"></param>
        /// <returns></returns>
        public double[] normpdf(List<double> datos, double mu, double sigma)
        {
            double exp;
            double fx;
            double[] y = new double[datos.Count];
            for (int i = 0; i < datos.Count(); i++)
            {
                 exp = ((-1)*Math.Pow((datos[i] - mu), 2)) / (2 * Math.Pow(sigma, 2));
                 fx = (1 / (sigma * Math.Sqrt(2 * Math.PI))) * Math.Exp(exp);
                 y[i] = fx;
            }
            return y;
        }
        /// <summary>
        /// Valor que te devuelve al evaluar la funcion en x 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="mu"></param>
        /// <param name="sigma"></param>
        /// <returns></returns>
        public double likelihood(double x, double mu, double sigma)
        {
            double exp = ((-1) * Math.Pow((x - mu), 2)) / (2 * Math.Pow(sigma, 2));
            double fx = (1 / (sigma * Math.Sqrt(2 * Math.PI))) * Math.Exp(exp);
            return fx;
        }
        /// <summary>
        /// Funcion que devuelve el valor de la integral por el método trapezoidal 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public double trapz(List<double> x, Func<double, double> y)
        {
            double r;
            int n = x.Count;
            double x1 = x[0];
            double x2 = x[x.Count()-1];
            double dx = (x2 - x1) / n;
            double s = 0.5 * (y(x1) + y(x2));
            for (int i = 1; i <= n; i++)
            {
                s = s + y(x1 + i*dx);
            }
            r = dx * s;
            return r;

        }
        /// <summary>
        /// devuele una tupla con los valores de m (pendiente de una recta) y b (la ordenada al origen) en ese orden de dos listas de datos ordenados 
        /// </summary>
        /// <param name="x"></param> lista 
        /// <param name="y"></param> lista
        /// <returns name = "m"></returns Tuple >
        /// <returns name = "b"></returns Tuple >
        public Tuple<double, double> regresion(List<double> x, List<double> y)
        {
            double m = 0, b = 0, prod_xy = 0, suma_x2 = 0;
            double n = x.Count;
            double suma_x = sum(x);
            double suma_y = sum(y);
            double x2 = Math.Pow(suma_x, 2);

            for (int i = 0; i < x.Count; i++)
            {
                prod_xy = prod_xy + (x[i] * y[i]);
            }
            for (int i = 0; i < x.Count; i++)
            {
                suma_x2 = suma_x2 + Math.Pow(x[i], 2);
            }

            m = ((n * prod_xy) - ((suma_x) * (suma_y))) / ((n * suma_x2) - x2);

            b = (suma_y - (m * suma_x)) / n;
            return new Tuple<double, double>(m, b);

        }
    }
}
