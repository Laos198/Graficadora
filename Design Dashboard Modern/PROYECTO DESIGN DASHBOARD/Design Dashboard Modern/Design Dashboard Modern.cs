using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace Design_Dashboard_Modern
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void Salir_Click(object sender, EventArgs e)
        {
            Application.Exit();
        }

        private void Minimizar_Click(object sender, EventArgs e)
        {
            WindowState = FormWindowState.Minimized;
        }

        private void Maximizar_Click(object sender, EventArgs e)
        {
            WindowState = FormWindowState.Maximized;
            Maximizar.Visible = false;
            Restaurar.Visible = true;
        }

        private void Restaurar_Click(object sender, EventArgs e)
        {
            WindowState = FormWindowState.Normal;
            Restaurar.Visible = false;
            Maximizar.Visible = true;
        }

        private void MenuSidebar_Click(object sender, EventArgs e)
        {
            if(Sidebar.Width == 270)
            {
                Sidebar.Visible = false;
                Sidebar.Width = 68;
                SidebarWrapper.Width = 90;
                LineaSidebar.Width = 52;
                AnimacionSidebar.Show(Sidebar);
            }
            else
            {
                Sidebar.Visible = false;
                Sidebar.Width = 270;
                SidebarWrapper.Width = 300;
                LineaSidebar.Width = 252;
                AnimacionSidebarBack.Show(Sidebar);
            }
        }

        private void Temporizador_Tick(object sender, EventArgs e)
        {
            Temporizador.Stop();

            this.Grafico.Series["ChartLinea"].Points.AddXY("LUNES", 10);
            this.Grafico.Series["ChartLinea"].Points.AddXY("MARTES", 20);
            this.Grafico.Series["ChartLinea"].Points.AddXY("MIERCOLES", 30);
            this.Grafico.Series["ChartLinea"].Points.AddXY("JUEVES", 40);
            this.Grafico.Series["ChartLinea"].Points.AddXY("VIERNES", 50);
            this.Grafico.Series["ChartLinea"].Points.AddXY("SABADO", 40);
            this.Grafico.Series["ChartLinea"].Points.AddXY("DOMINGO", 20);
        }
    }
}
