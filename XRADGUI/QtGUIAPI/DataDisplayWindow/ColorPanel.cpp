#include "pre.h"
#include "ColorPanel.h"

ColorPanel::ColorPanel(QWidget* parent)
	: QWidget(parent)
{
	m_main_layout = new QHBoxLayout(this);
	m_main_layout->setSpacing(5);
	m_main_layout->setMargin(0);

	m_title_label = new QLabel("Color:");
	m_red_label = new QLabel("R:");
	m_green_label = new QLabel("G:");
	m_blue_label = new QLabel("B:");


	m_red_spin_box = getNewSpinBox(0, 255);
	m_gren_spin_box = getNewSpinBox(0, 255);
	m_blue_spin_box = getNewSpinBox(0, 255);

	m_main_layout->addWidget(m_title_label);
	m_main_layout->addWidget(m_red_label);
	m_main_layout->addWidget(m_red_spin_box);
	m_main_layout->addWidget(m_green_label);
	m_main_layout->addWidget(m_gren_spin_box);
	m_main_layout->addWidget(m_blue_label);
	m_main_layout->addWidget(m_blue_spin_box);


	for (size_t i = 0; i < 5; i++)
	{
		m_items.push_back(new QRadioButton);
		m_colors.push_back(QColor(QColor::fromRgb(i % 2 == 0 && i != 4 ? 255 : 0, i % 3 == 0 ? 255 : 0, i % 4 == 0 ? 255 : 0)));
	}
	for (size_t i = 0; i < m_items.size(); i++)
	{

		m_items[i]->setStyleSheet(getStyleSheet(m_colors[i]));
		m_items[i]->setObjectName((QString)int(i));
		m_items[i]->setFocusPolicy(Qt::TabFocus);
		m_main_layout->addWidget(m_items[i]);
		connect(m_items[i], &QRadioButton::clicked, this, &ColorPanel::click);
	}



	m_items[0]->setChecked(true);
	setColorPickerValue(m_colors[0]);

	installEventFilter(this);
}
ColorPanel::~ColorPanel()
{

}

bool ColorPanel::eventFilter(QObject* target, QEvent* event)
{
	if (target == m_blue_spin_box || target == m_red_spin_box || target == m_gren_spin_box)
	{
   		setColorPanelValue(QColor::fromRgb(m_red_spin_box->value(), m_gren_spin_box->value(), m_blue_spin_box->value()));
	}
	return parent()->eventFilter(target, event);
}

QSpinBox* ColorPanel::getNewSpinBox(int min, int max)
{
	QSpinBox* result = new QSpinBox(this);
	result->setRange(min, max);
	result->setFocusPolicy(Qt::FocusPolicy::ClickFocus);

	result->installEventFilter(this);

	return result;
}
QColor ColorPanel::currentColor()
{
	return m_colors[m_changed_item_number];
}

void ColorPanel::setColorPanelValue(const QColor& in_color)
{
	m_colors[m_changed_item_number] = in_color;
	m_items[m_changed_item_number]->setStyleSheet(getStyleSheet(in_color));
}

void ColorPanel::setColorPickerValue(const QColor& color)
{
	m_red_spin_box->setValue(color.red());
	m_gren_spin_box->setValue(color.green());
	m_blue_spin_box->setValue(color.blue());
}

void ColorPanel::click()
{
	QObject* obj = sender();

	for (size_t i = 0; i < m_items.size(); i++)
	{
		if (obj->objectName() == m_items[i]->objectName())
		{
			m_changed_item_number = i;
			setColorPickerValue(m_colors[i]);
		}
	}

}

